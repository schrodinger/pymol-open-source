# -c

/print "BEGIN-LOG"

load dat/names.pdb

select test,name O4'
select test,*/O4'

select test,*/O4'+O3'
select test,(*/O4',O3')
select test,name O4'+O3'
select test,name "O4'+O3'"

select test,name O4'+Na\+
select test,(name Na\+,O4')
select test,name Na\++O4'
select test,*/Na\++O4'

select test,*/O4'+O4

select test,*/O2\*+O2

select test,*/O2\*+O2'

# ugly ugly ugly

select test,name 'O4''
select test,name "O4'"+O2'
select test,name 'O4''+O2'

select test,name O4'+C2+C4+Na\++Cl\-
select test,*/O4'+C2+C4+Na\++Cl\-

select test,names////O4'
select test,names and */O4'+'Na\+'+'C2'+'C4'

select test,chain ''
select test,segi ''

select test,chain R+''
select test,chain ''+R

select test,name Na\+
select test,name "Na\+"

/print "END-LOG"


