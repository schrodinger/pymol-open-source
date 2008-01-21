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

tmp_sele1 = "assign_tmp1"
tmp_sele2 = "assign_tmp2"

from chempy import champ
from chempy.champ import Champ
from pymol import cmd

def missing_c_termini(selection="(all)",quiet=0,_self=cmd):
    cmd=_self
    # assumes that hydogens are not present!
    
    sele_list = []
    ch=Champ()
    model = cmd.get_model(selection)
    model_pat = ch.insert_model(model)
    assn_pat = ch.insert_pattern_string("[N+0+1]C[C;D2]<0>(=O)")
    ch.pattern_clear_tags(model_pat)
    if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
        result = ch.pattern_get_ext_indices_with_tags(model_pat)
        for atom_tag in result[0]: # just iterate over atom tags
            if len(atom_tag[1])==1: # one and only one match
                if atom_tag[1][0]==0:
                    sele_list.append(atom_tag[0])
    cmd.select_list(tmp_sele1,selection,sele_list, mode='index')
    while cmd.pop(tmp_sele2,tmp_sele1)>0: # complete the carboxy terminus
        cmd.edit(tmp_sele2)
        cmd.attach("O",1,1,"OXT",quiet=1)
        cmd.unpick()
    cmd.delete(tmp_sele1)
    
    
def formal_charges(selection="(all)",quiet=0,_self=cmd):
    cmd=_self
    result = 1
    # assumes that hydogens are not present!
    
    # first, set all formal charges to zero
    
    cmd.alter(selection,"formal_charge=0")

    # next, flag all atoms so that we'll be able to detect what we miss
    
    cmd.flag(23,selection,'set')

    # get the residue dictionary for formal charges
    
    if not hasattr(champ,'formal_charge_dict'):
        from chempy.champ.formal_charges import formal_charge_dict
        champ.formal_charge_dict = formal_charge_dict

    # iterate through the residue dictionary matching each residue based on chemistry
    # and generating the expressions for reassigning formal charges
    
    alter_list = []
    for resn in champ.formal_charge_dict.keys():
        if cmd.select(tmp_sele1,"(%s) and resn %s"%(selection,resn))>0:
            entry = champ.formal_charge_dict[resn]
            for rule in entry:
                model = cmd.get_model(tmp_sele1)
                ch = Champ()
                model_pat = ch.insert_model(model)         
                assn_pat = ch.insert_pattern_string(rule[0])
                ch.pattern_clear_tags(model_pat)
                if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
                    result = ch.pattern_get_ext_indices_with_tags(model_pat)
                    for atom_tag in result[0]: # just iterate over atom tags
                        if len(atom_tag[1])==1: # one and only one match
                            tag = atom_tag[1][0]
                            formal_charge = rule[1][tag]
                            # the following expression both changes the formal charge and resets flag 23
                            alter_list.append([atom_tag[0],
                                                     "formal_charge=%d;flags=flags&0x-8388609"%formal_charge])

    if 1: # n-terminal amine
        # non-proline 
        ch=Champ()
        model = cmd.get_model(selection)
        model_pat = ch.insert_model(model)
        assn_pat = ch.insert_pattern_string("[N;D1]<0>CC(=O)")
        ch.pattern_clear_tags(model_pat)
        if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
                result = ch.pattern_get_ext_indices_with_tags(model_pat)
                for atom_tag in result[0]: # just iterate over atom tags
                    if len(atom_tag[1])==1: # one and only one match
                        if atom_tag[1][0]==0:
                            # the following expression both changes the formal charge and resets flag 23
                            alter_list.append([atom_tag[0],
                                                     "formal_charge=1;flags=flags&-8388609"])
        # proline residues
        ch=Champ()
        model = cmd.get_model(selection)
        model_pat = ch.insert_model(model)
        assn_pat = ch.insert_pattern_string("C1CC[N;D2]<0>C1C(=O)")
        ch.pattern_clear_tags(model_pat)
        if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
                result = ch.pattern_get_ext_indices_with_tags(model_pat)
                for atom_tag in result[0]: # just iterate over atom tags
                    if len(atom_tag[1])==1: # one and only one match
                        if atom_tag[1][0]==0:
                            # the following expression both changes the formal charge and resets flag 23
                            alter_list.append([atom_tag[0],
                                                     "formal_charge=1;flags=flags&-8388609"])
                                    
    if 1: # c-terminal acid
        ch=Champ()
        model = cmd.get_model(selection)
        model_pat = ch.insert_model(model)
        assn_pat = ch.insert_pattern_string("NCC(=O<0>)[O;D1]<1>")
        ch.pattern_clear_tags(model_pat)
        if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
                result = ch.pattern_get_ext_indices_with_tags(model_pat)
                for atom_tag in result[0]: # just iterate over atom tags
                    if len(atom_tag[1])==1: # one and only one match
                        if atom_tag[1][0]==1:
                            # the following expression both changes the formal charge and resets flag 23
                            alter_list.append([atom_tag[0],
                                                     "formal_charge=-1;flags=flags&-8388609"])
        
    # now evaluate all of these expressions efficiently en-masse 
    cmd.alter_list(selection,alter_list)

    # see if we missed any atoms
    missed_count = cmd.count_atoms("("+selection+") and flag 23")

    if missed_count>0:
        if not quiet:
            # looks like we did, so alter the user
            print " WARNING: %d atoms did not have formal charges assigned"%missed_count
        result = 0
    # remove the temporary selection we used to select appropriate residues
    
    cmd.delete(tmp_sele1)
    
    return result

def amber99(selection="(all)",quiet=0,_self=cmd):
    cmd=_self
    result = 1
    # first, set all parameters to zero

    cmd.alter(selection,"name=''")
    cmd.alter(selection,"partial_charge=0")
    cmd.alter(selection,"elec_radius=0.0")
    cmd.alter(selection,"text_type=''")

    # next, flag all atoms so that we'll be able to detect what we miss
    
    cmd.flag(23,selection,'set')

    # get the amber99 dictionary
    
    if not hasattr(champ,'amber99_dict'):
        from chempy.champ.amber99 import amber99_dict
        champ.amber99_dict = amber99_dict

    # iterate through the residue dictionary matching each residue based on chemistry
    # and generating the expressions for reassigning formal charges
    
    alter_list = []
    for resn in champ.amber99_dict.keys():
        if cmd.select(tmp_sele1,"(%s) and resn %s"%(selection,resn))>0:
            entry = champ.amber99_dict[resn]
            for rule in entry:
                model = cmd.get_model(tmp_sele1)
                ch = Champ()
                model_pat = ch.insert_model(model)
                ch.pattern_detect_chirality(model_pat)
                assn_pat = ch.insert_pattern_string(rule[0])
                ch.pattern_clear_tags(model_pat)
                if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
                    result = ch.pattern_get_ext_indices_with_tags(model_pat)
                    for atom_tag in result[0]: # just iterate over atom tags
                        if len(atom_tag[1])==1: # one and only one match
                            tag = atom_tag[1][0]
                            prop_list = rule[1][tag]
                            # the following expression both changes the formal charge and resets flag 23
                            alter_list.append([atom_tag[0],
            "name='''%s''';text_type='''%s''';partial_charge=%f;elec_radius=%f;flags=flags&-8388609"%prop_list])

    # now evaluate all of these expressions efficiently en-masse 
    cmd.alter_list(selection,alter_list)

    # see if we missed any atoms
    missed_count = cmd.count_atoms("("+selection+") and flag 23")

    if missed_count>0:
        if not quiet:
            # looks like we did, so alter the user
            print " WARNING: %d atoms did not have properties assigned"%missed_count
        result = 0

    # remove the temporary selection we used to select appropriate residues
    
    cmd.delete(tmp_sele1)

    return result

    
