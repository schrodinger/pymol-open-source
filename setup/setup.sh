#!/bin/bash
echo ' '
echo '============================================'
echo 'Creating "./pymol" startup script with '
echo "PYMOL_PATH=`pwd`"
echo '============================================'

touch pymol || exit 1

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
echo 'PYTHONPATH=$PYTHONHOME/lib/python2.7:$PYTHONPATH' >> ./pymol
echo 'PYTHONPATH=$PYTHONHOME/lib/python2.7/lib-tk:$PYTHONPATH' >> ./pymol
echo 'LD_LIBRARY_PATH=$PYMOL_PATH/ext/lib:$LD_LIBRARY_PATH' >> ./pymol
echo 'LANG=C' >> ./pymol
echo 'export LD_LIBRARY_PATH' >> ./pymol
echo 'export PYTHONHOME' >> ./pymol
echo 'export PYTHONPATH' >> ./pymol
echo 'export LANG' >> ./pymol
echo 'exec $PYMOL_PATH/pymol.exe "$@"' >> ./pymol
chmod 755 ./pymol

# remove extra libs from ext/lib
(cd ext/libextra && names=`ls` && cd ../lib && rm -f $names)

# add extra libs to ext/lib if there seems to be a problem like
# missing library or incompatible library
while true; do
    name=(`./pymol -cq 2>&1 | grep 'lib[^ /]*\.so' | sed 's!.*\(lib[^ /]*\.so[.0-9]*\).*!\1!'`)
    test -z "$name" && break
    if [ -e "ext/lib/$name" -o ! -e "ext/libextra/$name" ]; then
        echo "warning: problem with $name"
        break
    fi
    (cd ext/lib && ln -s ../libextra/$name .)
done

