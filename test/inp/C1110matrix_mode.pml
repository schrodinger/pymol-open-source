# -c

print "BEGIN-LOG"

reinit
set matrix_mode, 0

load dat/pept.pdb

remove not resi 5

iterate_state 1,all,print "%8.3f %8.3f %8.3f"%(x,y,z)

rotate x, 53, object=pept
rotate y, 27, object=pept
rotate z, 6, object=pept
translate [1,2,3], object=pept

save tmp/C1110a.pdb
for a in open("tmp/C1110a.pdb").readlines(): print a,

save tmp/C1110b.pdb, ref=pept
for a in open("tmp/C1110b.pdb").readlines(): print a,

load dat/pept.pdb, tmp
remove tmp and not resi 5

matrix_transfer pept, tmp

save tmp/C1110c.pdb, tmp,
for a in open("tmp/C1110c.pdb").readlines(): print a,

save tmp/C1110d.pdb, tmp, ref=tmp
for a in open("tmp/C1110d.pdb").readlines(): print a,

reinit
set matrix_mode, 1

load dat/pept.pdb

remove not resi 5

iterate_state 1,all,print "%8.3f %8.3f %8.3f"%(x,y,z)

rotate x, 53, object=pept
rotate y, 27, object=pept
rotate z, 6, object=pept
translate [1,2,3], object=pept

save tmp/C1110a.pdb
for a in open("tmp/C1110a.pdb").readlines(): print a,

save tmp/C1110b.pdb, ref=pept
for a in open("tmp/C1110b.pdb").readlines(): print a,

load dat/pept.pdb, tmp
remove tmp and not resi 5

matrix_transfer pept, tmp

save tmp/C1110c.pdb, tmp,
for a in open("tmp/C1110c.pdb").readlines(): print a,

save tmp/C1110d.pdb, tmp, ref=tmp
for a in open("tmp/C1110d.pdb").readlines(): print a,

reinit
set matrix_mode, 2

load dat/pept.pdb

remove not resi 5

iterate_state 1,all,print "%8.3f %8.3f %8.3f"%(x,y,z)

rotate x, 53, object=pept, object_mode=1
rotate y, 27, object=pept, object_mode=1
rotate z, 6, object=pept, object_mode=1
translate [1,2,3], object=pept, object_mode=1

save tmp/C1110a.pdb
for a in open("tmp/C1110a.pdb").readlines(): print a,

save tmp/C1110b.pdb, ref=pept
for a in open("tmp/C1110b.pdb").readlines(): print a,

load dat/pept.pdb, tmp
remove tmp and not resi 5

matrix_transfer pept, tmp

save tmp/C1110c.pdb, tmp,
for a in open("tmp/C1110c.pdb").readlines(): print a,

save tmp/C1110d.pdb, tmp, ref=tmp
for a in open("tmp/C1110d.pdb").readlines(): print a,


print "END-LOG"





