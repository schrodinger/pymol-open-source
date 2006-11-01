# -c

def dump44(lst): \
   new_lst = [] \
   for a in lst: \
      new_lst.append(int(a*1000)/1000.0) \
   print "[ %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[0:4]) \ 
   print "  %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[4:8]) \
   print "  %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[8:12]) \
   print "  %5.2f %5.2f %5.2f %5.2f ]"%tuple(new_lst[12:16])
 
/print "BEGIN-LOG"

load dat/pept.pdb,prot

map_new map1,gaussian,1.0,state=-3

translate [10,20,30], object=prot, object_mode=1
rotate x,30,prot,object=prot,object_mode=1
rotate y,30,prot,object=prot,object_mode=1
rotate z,30,prot,object=prot,object_mode=1

zoom prot

dump44(cmd.get_object_matrix("map1"))

matrix_transfer prot, map1

dump44(cmd.get_object_matrix("map1"))

isomesh m1, map1, 1.0, 5/, carve=2.1

ray renderer=2

isosurf s1, map1, 1.0, 5/, carve=2.1

ray renderer=2

translate [5,-10,2], object=prot, object_mode=1

rotate x,60,prot,object=prot,object_mode=1
rotate z,20,prot,object=prot,object_mode=1
rotate y,30,prot,object=prot,object_mode=1

matrix_transfer prot, map1

isomesh m1, map1, 1.0, 5/, carve=3

ray renderer=2

disable

isosurf s1, map1, 1.0, 5/, carve=3

ray renderer=2

disable

isomesh m2, map1, 1.0, 5/, 3

ray renderer=2

disable

isosurf s2, map1, 1.0, 5/, 3

ray renderer=2

disable

isomesh m3, map1

ray renderer=2

disable

isosurf s3, map1

ray renderer=2

enable

ray renderer=2

save tmp/C1200.pse
reinit
load tmp/C1200.pse
ray renderer=2

disable

isomesh m4, map1, 2.0, 2/, 4
isosurf s4, map1, 2.0, 2/, 4

ray renderer=2

disable

enable m1

isolevel m1, 3

ray renderer=2

disable

enable s1

isolevel s1, 5

ray renderer=2

print "END-LOG"