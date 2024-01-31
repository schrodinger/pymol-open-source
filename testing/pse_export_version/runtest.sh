#!/bin/bash
##
## Generate PyMOL Session Files (.pse) with different PyMOL versions
## and different "pse_export_version"s, and load them into previous
## program versions.
##

# exit on error
set -e

# where to find PyMOL installations
installprefix="/opt/pymol-"
installsuffix="*/pymol"

# temporary directory
tmpdir=.

# pse_binary_dump=0/1
bin_dump=1

gen() {
    testversions_outer=""
    while [[ $# > 0 ]]; do
        v_in=$1

        # pse_export_version new in PyMOL 1.7.6
        if [[ $v_in < 1.76 ]]; then
            return
        fi

        # for all installations of "v_in" version
        for pymol in $installprefix${v_in}${installsuffix}; do
            echo $pymol
            build=$(basename $(dirname $pymol))
            testversions="$testversions_outer"

            # for all applicable "pse_export_version"s
            for v_out in $*; do
                filename="$tmpdir/$build-ev$v_out.pse"

                # generate session file
                $pymol -ckq makesession.pml \
                    -d "set pse_binary_dump, $bin_dump" \
                    -d "set pse_export_version, $v_out" \
                    -d "save $filename"

                testversions="$v_out $testversions"

                # load session file with older installations
                for v_test in $testversions; do
                    for pymoltest in $installprefix${v_test}${installsuffix}; do
                        echo "TESTING $pymoltest $filename"
                        $pymoltest -ckq $filename verify.py
                    done
                done
            done
        done

        shift
        testversions_outer="$v_in $testversions_outer"
    done
}

# run on these versions
# (must be floats, and must match the installation directories)
gen 2.0 1.86 1.74 1.2

echo DONE
