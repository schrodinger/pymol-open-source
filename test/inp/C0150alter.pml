# -c

/print "BEGIN-LOG"

load dat/tiny.pdb

alter (all),resn = name
alter (all),b = b + 10
alter (all),q = 0.95
alter (all),name = name[0:2]
alter (het),type = 'ATOM'
alter (resi 9),type ='HETATM'
alter (e;c),elem = 'X'
alter (resi 8),resi=str(int(resi)+10)
alter (resi 7),alt='B'
save cmp/C0150alter.1.pdb

/print "END-LOG"
