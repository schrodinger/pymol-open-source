#!/bin/csh
/usr/bin/wish <<EOF >& $1
puts [ selection get ]
exit
EOF
if ( $status ) then
echo "" > $1
endif
