# pymol -x

load dat/3al1.pdb,t0
create t1 = (t0 & (alt '',A))
create t2 = (t1 & (i; 101:105))
create t3 = (t1 & (i; 104:111))
create t4 = (t2 | t3)
create t5 = (t0 & ((t0 & i;100) x;5)),1,1
create t5 = (t0 & ((t0 & i;100) x;7)),1,2
create t5 = (t0 & ((t0 & i;100) x;9)),1,3
create t5 = (t0 & ((t0 & i;100) x;11)),1,4
create t5 = (t0 & ((t0 & i;100) x;13)),1,5 
create t5 = (t0),1,3
save cmp/cmd_create01.01.pdb,(t1)
save cmp/cmd_create01.02.pdb,(t2)
save cmp/cmd_create01.03.pdb,(t3)
save cmp/cmd_create01.04.pdb,(t4)
save cmp/cmd_create01.05.pdb,(t5)
save cmp/cmd_create01.06.pdb,(t5),2
save cmp/cmd_create01.07.pdb,(t5),3
save cmp/cmd_create01.08.pdb,(t5),4
save cmp/cmd_create01.09.pdb,(t5),5
refresh
quit
