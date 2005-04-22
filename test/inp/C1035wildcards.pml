# -c

/print "BEGIN-LOG"

load dat/pept.pdb

count_atoms name C*
count_atoms */C*
count_atoms A*/
count_atoms A*/O*

alter all,segi='FOUR'
alter 1-5/,segi='FIRE'
alter 6-7/,segi='SORE'

count_atoms F*///
count_atoms F*R*///
count_atoms F*R///
count_atoms *R*///
count_atoms *IRE///
count_atoms *RE///
count_atoms *E///
count_atoms *O*///

count_atoms chain E

alter 1-3/,chain='A'
alter 6-9/,chain='B'
alter 10-11/,chain='G'

count_atoms chain A
count_atoms chain A:B
count_atoms chain B:F
count_atoms chain A:Z
count_atoms chain A+E
count_atoms chain A+F

count_atoms */*/*/*/*
count_atoms */*/*/*
count_atoms */*/*
count_atoms */*
count_atoms *

count_atoms /*/*/*/*/*
count_atoms /*/*/*/*
count_atoms /*/*/*
count_atoms /*/*
count_atoms /*

count_atom /*/*RE
count_atom /*////ca
count_atom *////ca
count_atom /////ca

create pap,pept and not elem o
create abt,pept and not elem n

count_atoms p*
count_atoms *t
count_atoms *pt
count_atoms p*p
count_atoms p*p*
count_atoms a*

count_atoms /p*/*RE
count_atoms /*t/*O*
count_atoms /*pt//A:E

count_atoms /p*/*RE/E/h*/c*

/print "END-LOG"









