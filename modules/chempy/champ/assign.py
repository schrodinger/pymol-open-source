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

tmp_sele = "_assing_tmp"

from chempy import champ
from chempy.champ import Champ
from pymol import cmd

def formal_charges(object):

   # first, set all formal charges to zero
   
   cmd.alter(object,"formal_charge=0")

   # next, flag all atoms so that we'll be able to detect what we miss
   
   cmd.flag(31,object,'set')

   # get the residue dictionary for formal charges
   
   if not hasattr(champ,'formal_charge_dict'):
      from chempy.champ.formal_charges import formal_charge_dict
      champ.formal_charge_dict = formal_charge_dict


   # iterate through the residue dictionary matching each residue based on chemistry
   # and generating the expressions for reassigning formal charges
   
   alter_list = []
   for resn in champ.formal_charge_dict.keys():
      if resn=='LEU':
         break
      if cmd.select(tmp_sele,"(%s) and resn %s"%(object,resn))>0:
         entry = champ.formal_charge_dict[resn]
         model = cmd.get_model(tmp_sele)
         ch = Champ()
         model_pat = ch.insert_model(model)         
         assn_pat = ch.insert_pattern_string(entry[0])
         ch.pattern_clear_tags(model_pat)
         if ch.match_1v1_n(assn_pat,model_pat,10000,2)>0:
            result = ch.pattern_get_ext_indices_with_tags(model_pat)
            for atom_tag in result[0]: # just iterate over atom tags
               if len(atom_tag[1])==1: # one and only one match
                  tag = atom_tag[1][0]
                  formal_charge = entry[1][tag]
                  # the following expression both changes the formal charge and resets flag 31
                  alter_list.append([atom_tag[0],
                                     "formal_charge=%d;flags=flags&0x7FFFFFFF"%formal_charge])

   # now evaluate all of these expressions efficiently en-masse 
   cmd.alter_list(object,alter_list)

   # see if we missed any atoms
   missed_count = cmd.count_atoms("flag 31")

   if missed_count>0:
      # looks like we did, so alter the user
      print "WARNING: %d atoms did not have formal charges assigned"%missed_count

   # remove the temporary selection we used to select appropriate residues
   
   cmd.delete(tmp_sele)
   

   
