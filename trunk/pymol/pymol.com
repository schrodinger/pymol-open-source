#!/bin/csh 
#
# PyMOL startup script
#
setenv PYMOL_PATH ${WLD_MTK}/pkg/pymol
setenv PYMOL_EXTLIBPATH $PYMOL_PATH/ext/lib
#
# Tcl/Tk path
setenv TCL_LIBRARY $PYMOL_EXTLIBPATH/tcl8.0
#
# dynamic linking
# 
if ( $?LD_LIBRARY_PATH ) then
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYMOL_EXTLIBPATH}
else
setenv LD_LIBRARY_PATH ${PYMOL_EXTLIBPATH}
endif
#
# python modules
#
if ( $?PYTHONPATH ) then
setenv PYTHONPATH ${PYMOL_EXTLIBPATH}/python1.5:${PYTHONPATH}
else
setenv PYTHONPATH ${PYMOL_EXTLIBPATH}/python1.5
endif
#
#
#dbx $PYMOL_PATH/ext/bin/python 
#gdb $PYMOL_PATH/ext/bin/python 
$PYMOL_PATH/ext/bin/python $PYMOL_PATH/modules/launch_pymol.py $*

