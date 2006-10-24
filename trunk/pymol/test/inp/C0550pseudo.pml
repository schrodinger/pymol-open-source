# -c

/print "BEGIN-LOG"

pseudoatom test
move x,10
pseudoatom test

count_atoms
ray renderer=2

iterate all,print name+ " %0.2f"%vdw

dele all

load dat/pept.pdb
pseudoatom test, pept and resi 1
pseudoatom test, pept and resi 2
pseudoatom test, pept and resi 10

iterate name PS*, print name+ " %0.2f"%vdw
hide
as spheres, test

ray renderer=2

bond name PS2, name PS3
show sticks, test
ray renderer=2

pseudoatom test, vdw=5.0
ray renderer=2
pseudoatom test, color=blue
ray renderer=2
pseudoatom test, label=Hello world
ray renderer=2


dele all

load dat/pept.pdb
pseudoatom pept, pept and resi 1
pseudoatom pept, pept and resi 2
pseudoatom pept, pept and resi 10

count_atoms pept

bond name PS2, name PS3

hide nonbonded
ray renderer=2

dele all

load dat/pept.pdb
pseudoatom test, all, b=1.5, q=1.6, resn=abc, resi=10, name=ilk, elem=s, vdw=5, segi=asdf, chain=''
iterate name ilk,print "%8.3f %8.3f %8.3f"%(b,q,vdw)

/print "END-LOG"

