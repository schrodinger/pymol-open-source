#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"
puts [ selection get ]
exit

