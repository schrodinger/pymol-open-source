# pymol -x

load dat/water.pdb,wat
load dat/helix_amber.pdb,prot

color blue,wat
color salmon,prot

select t1 = (byres wat and (prot gap -0.4))
create t2 = (byres wat and (prot gap -0.2))
create t3 = (byres wat and (prot gap  0.2))
create t4 = (byres wat and (prot gap  0.4))
create t5 = (byres wat and (prot gap  0.6))

color orange,(t1 &!(t1 in t4))
color yellow,(t1 &!(t1 in t3))
color red,(t1 &!(t1 in t2))
hide all
show lines,(prot | t1)
save cmp/cmd_gap01.01.pdb,(t3)
refresh
quit
