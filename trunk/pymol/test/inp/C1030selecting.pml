# -c

/print "BEGIN-LOG"

load dat/pept.pdb

select (resi trp)
indicate (name ca)
deselect
select 1/o
select 1/o+ca+n
select /pept

select bb = (name ca or name c or name n)
select bb = (name ca,c,n)
select bb = */ca+c+n
select bb = (*/ca,c,n)
select bb = (not (not (name ca or name c or name n)))
select bb = (not (not (n. c|n. ca|n. n)))
select bb = (not (not (not (not (not (not (not (not (not (not (not (not (not (not (*/ca,c,n)))))))))))))))
select bb = (all and all and all and (not (not all)) and (not none) and all and not (not (n;ca,c,n)))
select bb,(name ca or (n;c | (n. n) | (*/n)))
select bb, \
name ca \
or name c \
or name n

select bb = bb
select bb = not not bb
select cc = not bb
select bb = not cc

enable
enable bb
enable cc
disable cc
hide lines,cc
hide lines,bb
show (bb)
set ray_default_renderer=2
ray

/print "END-LOG"









