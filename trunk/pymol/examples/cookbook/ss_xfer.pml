
# ss_xfer.pml:  Copies secondary structure from a source protein
# over to a target protein.  Residue identifiers must match perfectly.

load 1t46.pdb, source

load 1pkg.pdb, target

# remove all secondary structure assignments

alter all, ss=''

# show unassigned cartoon

as cartoon

# assign sequence to source

dss source

# now read and store assignments based on residue #'s 

stored.ss = {}

iterate source and polymer and name ca, stored.ss[resi]=ss

# and assign them to the target protein

alter target and polymer and name ca, ss=stored.ss.get(resi,'')

# rebuild the cartoon

rebuild




