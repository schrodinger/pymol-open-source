# CSH Setup script for Unix 
#
setenv PYMOL_PATH $HOME/pymol
setenv PYMOL_EXTLIBPATH $HOME/pymol-ext/lib
if ( $?LD_LIBRARY_PATH ) then
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYMOL_EXTLIBPATH}
else
setenv LD_LIBRARY_PATH $PYMOL_EXTLIBPATH
endif
#
alias pymol $PYMOL_PATH/pymol

