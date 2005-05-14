# -c

def dump_mat(lst): \
   new_lst = [] \
   for a in lst: \
      new_lst.append(int(a*1000)/1000.0) \
   print "[ %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[0:4]) \ 
   print "  %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[4:8]) \
   print "  %5.2f %5.2f %5.2f %5.2f  "%tuple(new_lst[8:12]) \
   print "  %5.2f %5.2f %5.2f %5.2f ]"%tuple(new_lst[12:16])
   
def dump_lst(lst): \
   new_lst = [] \
   for a in lst: \
      new_lst.append(int(a*1000)/1000.0) \
   print "[", \
   for a in new_lst: \
     print "%5.2f"%a, \
   print "]" 
   
/print "BEGIN-LOG"

load dat/pept.pdb, mol1
load dat/pept.pdb, mol2
load dat/pept.pdb, mol3

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))

# move mol1

translate [2,3,4], object=mol2, object_mode=1

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))

# apply mol2's history matrix to mol3's coordinates

matrix= cmd.get_object_matrix("mol2")

cmd.transform_object("mol3",matrix, homogenous=1)

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))


# now try again using the matrix_transfer command

dele mol3
load dat/pept.pdb,mol3

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))

# but first, try the null case: mol1 hasn't been moved, so the following
# should have no effect

matrix_transfer mol1, mol3

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))


# whereas, this should move mol3 onto mol2

matrix_transfer mol2, mol3

dump_lst((cmd.get_distance("mol1///5/ca","mol2///5/ca"), \
 cmd.get_distance("mol1///5/ca","mol3///5/ca"), \
 cmd.get_distance("mol2///5/ca","mol3///5/ca")))

dump_mat(cmd.get_object_matrix("mol1"))
dump_mat(cmd.get_object_matrix("mol2"))
dump_mat(cmd.get_object_matrix("mol3"))


# okay, now confirm that history accumulation works

dele all
load ../data/tut/1hpv.pdb
create ligand, organic
create protein, polymer
delete 1hpv

dump_lst((cmd.get_distance("ligand////s1","/protein/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/protein/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /protein/1HPV/B/VAL`82/CB")))

translate [2,3,4], object=protein, object_mode=1
rotate x, 20, object=protein, object_mode=1
rotate y, 40, object=protein, object_mode=1
rotate z, 90, object=protein, object_mode=1
translate [4,5,6], object=protein, object_mode=1

dump_lst((cmd.get_distance("ligand////s1","/protein/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/protein/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /protein/1HPV/B/VAL`82/CB")))

matrix_transfer protein, ligand

dump_mat(cmd.get_object_matrix("protein"))
dump_mat(cmd.get_object_matrix("ligand"))

dump_lst((cmd.get_distance("ligand////s1","/protein/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/protein/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /protein/1HPV/B/VAL`82/CB")))

# see if we can use "fit" to get the same matrix

load ../data/tut/1hpv.pdb, copy
create copy_ligand, copy and organic
create copy_protein, copy and polymer
create copy_protein2, copy and polymer
delete copy

matrix_transfer protein, copy_protein2
fit copy_protein, protein

dump_mat(cmd.get_object_matrix("protein"))
dump_mat(cmd.get_object_matrix("copy_protein"))
dump_mat(cmd.get_object_matrix("copy_protein2"))

# now see if we can fit back to the original

load ../data/tut/1hpv.pdb, target

dump_lst((cmd.get_distance("ligand////s1","/target/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/target/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /target/1HPV/B/VAL`82/CB")))

# first try undoing the ligand transformation

dump_mat(cmd.get_object_matrix("ligand"))
matrix_transfer target, ligand
dump_mat(cmd.get_object_matrix("ligand"))

dump_lst((cmd.get_distance("ligand////s1","/target/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/target/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /target/1HPV/B/VAL`82/CB")))

# now see if we can fit the protein back onto the original coordinates for a similar result

dump_mat(cmd.get_object_matrix("protein"))
fit protein, target

dump_mat(cmd.get_object_matrix("protein"))
dump_mat(cmd.get_object_matrix("ligand"))

dump_lst((cmd.get_distance("ligand////s1","/protein/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/protein/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /protein/1HPV/B/VAL`82/CB")))

# since we're back to the original coordinates, then something like
# this should have NO effect

matrix_transfer protein, ligand

dump_lst((cmd.get_distance("ligand////s1","/protein/1HPV/B/ILE`84/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/C4","/protein/1HPV/B/ILE`50/CB"), \
 cmd.get_distance("/ligand/1HPV//478`200/O3"," /protein/1HPV/B/VAL`82/CB")))

# try reset_matrix

dump_lst((cmd.get_distance("/target/1HPV/B/ILE`84/CB","/copy_protein2/1HPV/B/ILE`84/CB"),))
dump_mat(cmd.get_object_matrix("copy_protein2"))

matrix_reset copy_protein2
dump_mat(cmd.get_object_matrix("copy_protein2"))

dump_lst((cmd.get_distance("/target/1HPV/B/ILE`84/CB","/copy_protein2/1HPV/B/ILE`84/CB"),))

/print "END-LOG"