# This example mutates every residue of a chain to alanine.

# load a structure

cmd.load("$TUT/1hpv.pdb")

# activate the mutagensis wizard

cmd.wizard("mutagenesis")

# set target residue type

cmd.get_wizard().set_mode("ALA")

from pymol import stored
stored.list = []

# generate our list of residues (CA-atoms) to mutate

cmd.iterate("1hpv//A//CA",
"stored.list.append('/'.join([model,segi,chain,resi,name]))")

# now iterate through each residue

for sele in stored.list:
   
   # pick a residue

   cmd.edit(sele)

   # notify the wizard about the picked residue

   cmd.get_wizard().do_pick(0)

   # apply the mutation

   cmd.get_wizard().apply()

   # update the screen, showing the mutation

   cmd.indicate("byres "+sele)
   cmd.refresh()
   
# now close the wizard

cmd.set_wizard()
cmd.deselect()
