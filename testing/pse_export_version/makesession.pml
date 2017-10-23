wizard mutagenesis
bg yellow
load ~/fetch_path/1rx1.pdb
load ~/fetch_path/1rx1_2fofc.ccp4
as cartoon
color forest, polymer
color red, organic
color orange, inorganic
show sticks, organic
show spheres, inorganic
isomesh mesh, 1rx1_2fofc
color blue, mesh
scene F1, store, message=First Scene
