# -c

/print "BEGIN-LOG"

load dat/pept.pdb
load dat/pept.pdb
load dat/pept.pdb

get line_width
get line_width,doesnotexist
get line_width,doesnotexist,5
get line_width,pept
get line_width,pept,0
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3
get line_width,pept,4

set line_width,2.0,pept
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

set line_width,4.0,pept,1
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

set line_width,6.0,pept,2
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

unset line_width,pept,2
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

unset line_width,pept
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

unset line_width,pept,1
get line_width,pept
get line_width,pept,1
get line_width,pept,2
get line_width,pept,3

/print "END-LOG"

