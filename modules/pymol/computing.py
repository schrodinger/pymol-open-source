#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

import sys
cmd_module = __import__("sys").modules["pymol.cmd"]

import pymol

def clean(selection, present='', state=-1, fix='', restrain='',
          method='mmff', async_=0, save_undo=1, message=None,
          _self=cmd_module, **kwargs):
    '''
DESCRIPTION

    Run energy minimization on the given selection, using an MMFF94
    force field.
    '''
    raise pymol.IncentiveOnlyException
