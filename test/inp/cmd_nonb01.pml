# pymol -x

load dat/3al1.pdb,prot

show nonbonded
show nb_spheres

dist (id 597),(bonded),3.2

create t2 = (!bonded)
create t3 = (bonded)

refresh
quit
