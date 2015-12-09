'''
Various morphing workflows

(c) 2013 Schrodinger, Inc.
(c) 2011 Thomas Holder
(c) 2009 DeLano Scientific LLC.
'''

import pymol
cmd = __import__("sys").modules["pymol.cmd"]

def morph(name, sele1, sele2=None, state1=-1, state2=-1, refinement=3,
        steps=30, method='rigimol', match='align', quiet=1, _self=cmd):
    '''
DESCRIPTION

    Creates an interpolated trajectory between two or multiple conformations.
    If the two input objects are not the same, match them based on sequence
    alignment.

    This command supports two methods: rigimol and linear. RigiMOL is an
    incentive feature and only available to official PyMOL sponsors. Linear
    morphing is quick and robust but likely to produce distorted intermediates.

ARGUMENTS

    name = string: name of object to create

    sele1 = string: atom selection of first conformation

    sele2 = string: atom selection of second conformation {default: <sele1>}

    state1 = int: sele1 state {default: 1}. If state1=0 and sele1 has N
    states, create N morphings between all consecutive states and back from
    state N to 1 (so the morph will have N*steps states). If state2=0, create
    N-1 morphings and stop at last state.

    state2 = int: sele2 state {default: 2 if sele1=sele2, else 1}

    refinement = int: number of sculpting refinement cycles to clean
    distorted intermediates {default: 3}

    steps = int: number of states for sele2 object {default: 30}

    method = string: rigimol or linear {default: rigimol}

EXAMPLE

    fetch 1akeA 4akeA, async=0
    align 1akeA, 4akeA
    morph mout, 1akeA, 4akeA
    '''
    raise pymol.IncentiveOnlyException()

# vi: ts=4:sw=4:smarttab:expandtab
