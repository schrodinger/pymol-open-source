# -c

set auto_number_selections
/print "BEGIN-LOG"

load dat/tiny.pdb

load dat/il2.pdb

print cmd.get_names()

select tst = (all)
select tst = (none)
select tst = (name c)
dele tst
select tst = (segi e)
select tst = (all)

select (segi '')
select (name ca,c,n,o)
select (n;ca,c)
select (resn pro)
select (r;pro)
select (chain c,e)
select (c;c,e)
select (resi 100)
select (i;100)
select (id 50)
select (byres id 50)
select (alt '')

select (tiny)
select (tiny or il2)
select (tiny | il2)

select (vis)
hide (tiny)
select (vis)

print cmd.get_names('all')
delete sel*
print cmd.get_names('all')

select foo = (name c)
select fie = (name n,c)
select (foo and not fie)
select (fie &! foo)

select foe = (tiny a; 30)
select fum = (il2 and not foe)

select (hydro)
select (b<30 and b>0.1)
select (q=0.0)
select (q=1.0)

select (b=0.0)

load dat/pept.pdb

select tst,pept within 5 of tiny
select tst,pept within 50 of tiny

select tst, pept w. 1 of pept

print cmd.select("(none)")

dele all
load dat/pept.pdb,mult
load dat/il2.pdb,mult

frame 1
select tst,mult and present
frame 2
select tst,mult and present
select tst,state 1
select tst,state 2
hide
show sph,name c
select tst,present and vis
frame 1
select tst,present and vis


dele all
load dat/pept.pdb
indicate name ca
count_atoms indicate
indicate name *c*
count_atoms indicate
indicate (name *c*,*n*)
count_atoms indicate
indicate (name *c*+*n*)
count_atoms indicate
indicate name *G
count_atoms indicate
indicate */*G
count_atoms indicate

indicate rank -3+5+10-23+7+40-
count_atoms indicate

indicate id -3+5+10-23+7+40-
count_atoms indicate

indicate index -3+5+10-23+7+40-
count_atoms indicate

dele all
load dat/1tii.pdb
indicate A:C+F/10+90+100-/ca+c+n*

count_atoms indicate

count_atoms 1TII
count_atoms 1tii


count_atoms HOH`307/ expand 5

# implicit OR

count_atoms a// c//
count_atoms a// c// d//
count_atoms (a// c//) d//
count_atoms a// c//&none
count_atoms c//&none a//
count_atoms c//&none (a//&none a//)
count_atoms (c//&none) a// (a//&none)
count_atoms (a// (d// c//)) & (e// (b// a//))

/print "END-LOG"
