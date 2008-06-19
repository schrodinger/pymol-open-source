# This script illustrates how you can use PyMOL to count the
# number of molecules in a selection

cmd.load("$PYMOL_PATH/test/dat/1tii.pdb")

input_sele = "1tii and polymer"

prefix = "tmp_cnt_"

iter_sele = prefix + "iter"
mole_sele = prefix + "mole"

mole_cnt = 0

if cmd.select(iter_sele, input_sele)>0:
    while 1:
        atom_cnt = cmd.select(mole_sele, "bymol (first "+iter_sele+")")
        togo_cnt = cmd.select(iter_sele, iter_sele+" and not "+mole_sele)
        if atom_cnt>0:
            mole_cnt = mole_cnt + 1
            print "molecule %d contains %d atoms"%(mole_cnt,atom_cnt)
        if togo_cnt<=0:
            break
        

print "(%s) contains %d molecules total."%(input_sele, mole_cnt)

