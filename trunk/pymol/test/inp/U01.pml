# -c

# unix-specific externing tests

# Tests the commands:
#
# cd
# ls
# pwd
# system

set stop_on_exceptions=0

/print "BEGIN-LOG"

cd /
pwd
cd $PYMOL_PATH/test
ls
ls *.py
ls no_files
ls *no_files
system
system echo "hello"
system ls

/print "END-LOG"





