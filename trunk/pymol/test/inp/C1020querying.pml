# -c

#     count_atoms,        
#     count_states,       
#     dist,               
#     distance,           
#     export_dots,        
#     find_pairs,         
#     get_area,           
#     get_color_indices,  
#     get_color_tuple,    
#     get_dihedral,       
#     get_extent,         
#     get_model,          
#     get_names,          
#     get_names_of_type,  
#     get_phipsi,         
#     get_position,       
#     get_type,           
#     id_atom,            
#     identify,           
#     index,              
#     overlap,            
#     phi_psi

/print "BEGIN-LOG"

dele all
load dat/pept.pdb
count_atoms
count_atoms name ca
count_states
create pept,(all),1,2
create pept,(all),1,3
create pept,(none),1,2
create cpy,(all),1,1
count_states
count_states pept
count_states cpy
count_frames
mset 1 x10
count_frames

# these shouldn't produce output...

cmd.count_frames()
cmd.count_states()
cmd.count_atoms()

# distance

dist
distance
select lb,resi 5 and name ca
dist
select rb,resi 8 and name ca
dist
dist (lb),resi 10 and name ca
distance resi 6 and name n,resi 9 
distance resi 7 and name ca,all,3.0,zoom=1
print "%6.3f"%cmd.distance("resi 1 and name ca","resi 13 and name ca")
print "%6.3f"%cmd.distance("resi 1 and name ca","resi 14 and name ca")

# find_pairs

print cmd.find_pairs("name o","name n",mode=1,cutoff=3.4)

# get_area

dele all
load dat/pept.pdb 
get_area (none)
get_area
get_area state=2
get_area name ca
load dat/3al1.pdb
get_area
dele all

# get_color_indices

/lst=cmd.get_color_indices()[0:10]
print lst

# get_color_tuple

print cmd.get_color_tuple("red")
print cmd.get_color_tuple("blue")
print cmd.get_color_tuple("green")
print cmd.get_color_tuple("white")

#

dele all
load dat/pept.pdb
get_dihedral 1/n,1/ca,1/c,1/o

/lst=cmd.get_extent()
/cmd._dump_floats(lst[0]+lst[1])
/lst=cmd.get_extent("name ca")
/cmd._dump_floats(lst[0]+lst[1])

print len(cmd.get_model().atom)
print len(cmd.get_model("resi 1").bond)

# 

dele all
print cmd.get_names()
load dat/pept.pdb
print cmd.get_names()
edit 4/ca
print cmd.get_names("all")
print cmd.get_names("selections")
print cmd.get_names("objects")

print cmd.get_names_of_type("object:molecule")
print cmd.get_names_of_type("object:map")

# 
/pp=cmd.get_phipsi("resi 2")
kee = pp.keys()[0]
print kee

cmd._dump_floats(pp[kee])

# 

reset
get_position
zoom resi 10
get_position

#

get_type
get_type nonexistent
get_type pept

#

id_atom none
id_atom all
id_atom 5/ca

#

identify name ca

# 

index name ca

#

overlap pept,pept

alter (all),vdw=0.0
overlap pept,pept

alter (all),vdw=0.5
overlap pept,pept

#

phi_psi all

# 

dele all

load dat/pept.pdb
get_symmetry
get_symmetry pept
get_symmetry none
load dat/3al1.pdb
get_symmetry
get_symmetry 3al1
get_symmetry pept

hide
show cell,pept
ray renderer=2

apply(cmd.set_symmetry,["pept"]+cmd.get_symmetry("3al1"))

get_symmetry pept
get_symmetry 3al1
show cell,pept
ray renderer=2

/print "END-LOG"















