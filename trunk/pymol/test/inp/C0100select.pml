# -c

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

/print "END-LOG"
