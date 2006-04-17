
# load structures

load 1btu.pdb, 1btu
load 1ce5.pdb, 1ce5
load 1vgc.pdb, 1vgc

# separate the ligands

extract 1btu_lig, 1btu and organic
extract 1ce5_lig, 1ce5 and organic
extract 1vgc_lig, 1vgc and organic

# perform the alignment in the 1btu frame of reference

align 1ce5////ca, 1btu////ca
align 1vgc////ca, 1btu////ca

# bring the ligands along too 

matrix_transfer 1ce5, 1ce5_lig
matrix_transfer 1vgc, 1vgc_lig

# save out the coordinates in the 1ce5 frame of reference

save 1ce5_1ce5.pdb, 1ce5, ref=1ce5
save 1ce5_1btu.pdb, 1btu, ref=1ce5
save 1ce5_1vgc.pdb, 1vgc, ref=1ce5

save 1ce5_1ce5_lig.pdb, 1ce5_lig, ref=1ce5
save 1ce5_1btu_lig.pdb, 1btu_lig, ref=1ce5
save 1ce5_1vgc_lig.pdb, 1vgc_lig, ref=1ce5

