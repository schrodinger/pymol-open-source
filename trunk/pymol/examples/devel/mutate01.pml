# This example shows how to make a single mutation using a script
# note that no conformational analysis is done to insure
# that the resulting pose is at all reasonable.

# load a structure

load $PYMOL_PATH/test/dat/pept.pdb

# activate the mutagensis wizard

wizard mutagenesis

# set target residue type

cmd.get_wizard().set_mode("ILE")

# pick a residue

cmd.edit("pept///10/ca")

# notify the wizard about the picked residue
# (this will generate the mutation object)

cmd.get_wizard().do_pick(0)

# apply the mutation
# (this will delete the old residue, 
# substitue the new one, and form
# the peptide bonds) 

cmd.get_wizard().apply()

# now close the wizard

cmd.set_wizard()
