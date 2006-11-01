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

def attach_fragment(selection,fragment,hydrogen,anchor):
    if not selection in cmd.get_names("selections"):
        if fragment in cmd.get_names("objects"):
            print " Error: an object with than name already exists"
            raise QuietException
        else:
            cmd.fragment(fragment)
            if cmd.get_setting_legacy("auto_remove_hydrogens"):
                cmd.remove("(hydro and %s)"%fragment)
    else:
        cmd.fragment(fragment,tmp_editor)
        if cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
            cmd.fuse("(%s and id %d)"%(tmp_editor,hydrogen),"(pk1)",1)
            if cmd.get_setting_legacy("auto_remove_hydrogens"):
                cmd.remove("(hydro and pkmol)")            
        else:
            cmd.remove("(%s and id %d)"%(tmp_editor,hydrogen))
            cmd.fuse("(%s and id %d)"%(tmp_editor,anchor),"(pk1)",1)
            if cmd.get_setting_legacy("auto_remove_hydrogens"):
                cmd.remove("(hydro and pkmol)")            
        cmd.delete(tmp_editor)
        
def attach_amino_acid(selection,amino_acid,center=0,animate=-1):
    if (selection == 'pk1') and (selection not in cmd.get_names('all')):
        # create new object 
        if amino_acid in cmd.get_names("objects"):
            print " Error: an object with than name already exists"
            raise QuietException
        cmd.fragment(amino_acid)
        if cmd.get_setting_legacy("auto_remove_hydrogens"):
            cmd.remove("(hydro and %s)"%amino_acid)
        if cmd.count_atoms("((%s) and name c)"%amino_acid,quiet=1):
            cmd.edit("((%s) and name c)"%amino_acid)
        elif cmd.count_atoms("((%s) and name n)"%amino_acid,quiet=1):
            cmd.edit("((%s) and name n)"%amino_acid)            
    else:
        ss = int(cmd.get_setting_legacy("secondary_structure"))
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
        cmd.fragment(amino_acid,tmp_editor)
        if cmd.count_atoms("((%s) and elem n)"%selection,quiet=1):
            cmd.select(tmp_ed_save,"(%s)"%selection)
            cmd.iterate("(%s)"%selection,"stored.resv=resv")
            pymol.stored.resi = str(pymol.stored.resv-1)
            cmd.alter(tmp_editor,"resi=stored.resi")
            cmd.fuse("(%s and name C)"%(tmp_editor),tmp_ed_save,2)
            if cmd.get_setting_legacy("auto_remove_hydrogens"):
                cmd.remove("(pkmol and hydro)")
            if ((cmd.count_atoms("(name ca,ch3 and neighbor ?pk1)")==1) and
                (cmd.count_atoms("(name ca and neighbor ?pk2)"))):
                cmd.set_dihedral("(name ca and neighbor pk2)",
                                 "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
            if (cmd.select(tmp1,"?pk1")==1) and (cmd.select(tmp2,"?pk2")==1):
            
                if 0:
                    if ((cmd.select(tmp3,"(name ca and neighbor "+tmp2+"))")>0) and
                        (cmd.select(tmp4,"(name ca and neighbor "+tmp1+")")>0)):
                        cmd.set_dihedral( # PHI
                            tmp3, # CA +0
                            tmp2, # N +0 
                            tmp1, # C -1
                            tmp4, # CA -1
                            180.0) # insure that the peptide is planer

                cmd.set_geometry(tmp2,3,3) # make nitrogen planer
    
                if ss:
                    if amino_acid[0:3]!='pro':
                        if ((cmd.select(tmp4,
                                        "((!r;pro) and name c and neighbor (name ca and neighbor "
                                        +tmp2+"))")==1) and
                            (cmd.select(tmp3,"((!r;pro) and name ca and neighbor "+tmp2+")")==1)):
                            cmd.set_dihedral( # PHI
                                tmp4, # C
                                tmp3, # CA 
                                tmp2, # N
                                tmp1, # C
                                phi)

                    if ((cmd.select(tmp4,"(name n and neighbor (name ca and neighbor "+tmp1+"))")==1) and
                        (cmd.select(tmp3,"(name ca and neighbor "+tmp1+")")==1)):
                        cmd.set_dihedral( # PSI (n-1)
                            tmp2, # N
                            tmp1, # C
                            tmp3, # CA
                            tmp4, # N
                            psi)
                cmd.h_fix(tmp2) # fix hydrogen position
                        
            cmd.delete(tmp1)
            cmd.delete(tmp2)
            cmd.delete(tmp3)
            cmd.delete(tmp4)
                
            sele = ("(name N and (byres neighbor %s) and not (byres %s))"%
                      (tmp_ed_save,tmp_ed_save))
            if cmd.count_atoms(sele,quiet=1):
                cmd.edit(sele)
                cmd.center(sele,animate=animate)
            cmd.delete(tmp_ed_save)
                    
        elif cmd.count_atoms("((%s) and elem c)"%selection,quiet=1):
            cmd.select(tmp_ed_save,"(%s)"%selection)
            cmd.iterate("(%s)"%selection,"stored.resv=resv")
            pymol.stored.resi = str(pymol.stored.resv+1)
            cmd.alter(tmp_editor,"resi=stored.resi")
            cmd.fuse("(%s and name N)"%(tmp_editor),tmp_ed_save,2)
            if cmd.get_setting_legacy("auto_remove_hydrogens"):
                cmd.remove("(pkmol and hydro)")
            if (cmd.count_atoms("(name ca,ch3 and neighbor ?pk1)") and
                cmd.count_atoms("(name ca,ch3 and neighbor ?pk2)")):
                cmd.set_dihedral("(name ca,ch3 and neighbor pk2)",
                                 "(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
            if "pk1" in cmd.get_names('selections'):
                cmd.set_geometry("pk1",3,3) # make nitrogen planer
                cmd.h_fix("pk1") # fix hydrogen position
            if ss:
                if (cmd.select(tmp1,"?pk1")==1) and (cmd.select(tmp2,"?pk2")==1):
                    if amino_acid[0:3]=='nhh': # fix amide hydrogens
                        print "here1"
                        if ((cmd.select(tmp3,"(name h1 and neighbor "+tmp1+")")==1) and
                            (cmd.select(tmp4,"(name o and neighbor "+tmp2+")")==1)):
                            print "here2"
                            cmd.set_dihedral(
                                tmp4, # O
                                tmp2, # C
                                tmp1, # N
                                tmp3, # H1
                                180)
                    if amino_acid[0:3]!='pro':
                        if ((cmd.select(tmp3,"(name ca and neighbor "+tmp1+")")==1) and
                            (cmd.select(tmp4,"(name c and neighbor (name ca and neighbor "+tmp1+"))")==1)):
                            cmd.set_dihedral( # PHI
                                tmp2, # C
                                tmp1, # N
                                tmp3, # CA 
                                tmp4, # C
                                phi)
                    if ((cmd.select(tmp3,"(name ca and neighbor "+tmp2+")")==1) and
                        (cmd.select(tmp4,"(name n and neighbor (name ca and neighbor "+tmp2+"))")==1)):
                        cmd.set_dihedral( # PSI (n-1)
                            tmp4, # N
                            tmp3, # CA
                            tmp2, # C
                            tmp1, # N
                            psi)
                cmd.delete(tmp1)
                cmd.delete(tmp2)
                cmd.delete(tmp3)
                cmd.delete(tmp4)                               
            sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                      (tmp_ed_save,tmp_ed_save))
            if cmd.count_atoms(sele,quiet=1):
                cmd.edit(sele)
                cmd.center(sele,animate=animate)                
            cmd.delete(tmp_ed_save)
        elif cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
            print " Error: please pick a nitrogen or carbonyl carbon to grow from."
            cmd.delete(tmp_editor)
            raise QuietException
        else:
            print " Error: unable to attach fragment."
    cmd.delete(tmp_editor)
        

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
        
