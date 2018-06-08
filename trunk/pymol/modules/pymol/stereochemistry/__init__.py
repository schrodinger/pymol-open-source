import sys
import os
import pymol
cmd = sys.modules["pymol.cmd"]


def assign_stereo(selection='all', state=-1, method='', quiet=1, prop='stereo',
        _self=cmd):
    '''
DESCRIPTION

    Assign "stereo" atom property (R/S stereochemistry).

    Requires either a Schrodinger Suite installation (SCHRODINGER
    environment variable set) or RDKit (rdkit Python module).

USAGE

    assign_stereo [selection [, state [, method ]]]

ARGUMENTS

    selection = str: atom selection {default: all}

    state = int: object state {default: -1 (current)}

    method = schrodinger or rdkit: {default: try both}
    '''
    raise pymol.IncentiveOnlyException()
