# CSH Setup script for Unix 
#
setenv PYMOL_PATH $HOME/pymol
setenv PYMOL_EXTLIBPATH $HOME/pymol-ext/lib
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
alias pymol $PYMOL_PATH/pymol

