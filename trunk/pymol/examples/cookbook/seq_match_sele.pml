# this script uses the alignment engine to match aligned C-alpha atoms
# in related structures and graft a selection from one to the other

# load structures

fetch 1t46 1oky, async=0
as ribbon

# specify mobile and target objects / atoms

select mobile, 1t46////CA

select target, 1oky////CA

# specify the input selection to match across the alignment

select match_inp, 1oky///86-91/CA

##### generic from here down #####

# initialize the output selection

select match_out, none

# perform the alignment

align mobile, target, object=aln

python

# work around indexing bug in pre-1.2 versions

if cmd.get_version()[1]<1.2:
    fixver = 1
else:
    fixver = 0

# get raw alignment from Python

aln_list=cmd.get_raw_alignment("aln")

# iterate through matching atom sets

for match in aln_list:

    # init. temporary selections for accumulating matched atoms

    cmd.select("tmp_mob","none")
    cmd.select("tmp_tgt","none")

    # find corresponding atoms in mobile & target selectionsw

    for atom in match:
        atom_sele = "%s`%d"%(atom[0], atom[1]+fixver)
        cmd.select("tmp_mob","tmp_mob or (mobile and %s)"%atom_sele)
        cmd.select("tmp_tgt","tmp_tgt or (target and %s)"%atom_sele)

    # when both are found, match from input selection to output selection

    if (cmd.count_atoms("tmp_mob")==1) and (cmd.count_atoms("tmp_tgt")==1):
        if cmd.count_atoms("tmp_mob and match_inp"):
            cmd.select("match_out","match_out or tmp_tgt")
        elif cmd.count_atoms("tmp_tgt and match_inp"):
            cmd.select("match_out","match_out or tmp_mob")

python end

deselect
enable match_out
orient match_out
zoom match_out, 3