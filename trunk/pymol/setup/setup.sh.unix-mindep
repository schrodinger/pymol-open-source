#!/bin/sh
echo ' '
echo '============================================'
echo 'Creating "./pymol" startup script with '
echo "PYMOL_PATH=`pwd`"
echo '============================================'
echo 'If you need to move PyMOL in the future,'
echo 're-run "./setup.sh" afterwards.'
echo '============================================'
echo 'You may now want to copy or link:'
echo "\"`pwd`/pymol\"" 
echo 'to something like "/usr/local/bin/pymol"'
echo '============================================'
echo 'Enjoy!'
echo ' '
echo '#!/bin/sh' > pymol
echo '#' >> pymol
echo '# PyMOL startup script' >> ./pymol
echo '#' >> ./pymol
echo '# Set PYMOL_PATH to point to this directory' >> ./pymol
echo '#' >> ./pymol
echo '# ================================================' >> ./pymol
echo "PYMOL_PATH=`pwd`" >> ./pymol
echo '# ================================================' >> ./pymol
echo '#' >> ./pymol
echo 'export PYMOL_PATH' >> ./pymol
echo 'PYTHONHOME=$PYMOL_PATH/ext' >> ./pymol
echo 'export PYTHONHOME' >> ./pymol
echo 'exec $PYMOL_PATH/pymol.exe "$@"' >> ./pymol
chmod 755 ./pymol



