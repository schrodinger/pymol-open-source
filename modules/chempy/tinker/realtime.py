#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Scott Dixon, Metaphorics, LLC
#-*
#-*
#Z* -------------------------------------------------------------------

from __future__ import print_function

# pymol

from pymol import cmd

from chempy import io
from chempy import protein,hetatm
from chempy import feedback
from chempy import protein_amber99
from chempy import tinker
from chempy.tinker import keyword
from chempy.tinker.amber import Parameters,Topology,Subset
from chempy.tinker.state import State

import os

state = None
model = None

def assign(sele,preserve=0):

    from molobj import MolObj
    from typer import Typer,Rules
    import rules

    global state
    global model

    result = 1

    state = State()

    model = cmd.get_model(sele)
    # now assign atom types

    ruleSet = Rules()

#   ruleSet.fromList(rules.amber_types)
    ruleSet.fromList(rules.simple_types)

    mobj = MolObj()
    mobj.fromChemPyModel(model)

    typed = Typer(molObj = mobj)

    print(" realtime: assigning atom types")
    typed.applyRules(ruleSet)

    c = 0
    for a in typed.getNamedTypes():
        at = model.atom[c]
        if (at.text_type == '??') or (not preserve):
            if a=='':
                print(" warning: unable to assign atom type to atom %d"%c)
                result = 0
            else:
                cmd.alter("((%s) and (index %s))" % (sele,at.index),
                             "text_type ='%s'" % a)
                if feedback['tinker']:
                    print(" "+str(__name__)+': %s is a %s' % (at.name,a))
                at.text_type = a
        c = c + 1

    sm = 0
    for a in model.atom:
        a.resi = str(a.resi_number)
        sm = sm + a.partial_charge

    print(" lig: net charge on ligand  is %8.4f\n" % sm)

    return result

#   param = Parameters(tinker.params_path+"parm99_wld.dat")
#   param = Parameters(tinker.params_path+"simple_parm.dat")
#   param = Parameters("simple_parm.dat")

def setup(sele,preserve=0):

    global state
    global model

    state = State()

    model = cmd.get_model(sele)

    sm = 0
    for a in model.atom:
        a.resi = str(a.resi_number)
        sm = sm + a.partial_charge

    print(" lig: net charge on ligand  is %8.4f\n" % sm)

#   param = Parameters(tinker.params_path+"parm99_wld.dat")
    param = Parameters(tinker.params_path+"parm99_simple.dat")
    topo = Topology(model)

    subset = Subset(param,topo)

    if(subset.complete()):
        subset.write_tinker_prm("realtime.prm")

        state.params = "realtime.prm"

        state.load_model(model)
        return 1
    else:
        subset.dump_missing()
        model = model
        state = None
        return 0

def dyna(steps,iter=1):

    global state

    if not state:
        print(" realtime.dyna: please run setup first.")
    else:
        state.echo = 0

        model = state.model

        print(" realtime.dyna: %d atoms total\n" % model.nAtom)

        xtra_kw = []

        xtra_kw.extend(keyword.get_inactive(model,3))
        xtra_kw.extend(keyword.get_restrain_positions(model,2,0,10))
        xtra_kw.extend(keyword.get_restrain_positions(model,5,0.5,1))
        xtra_kw.extend(keyword.get_inactive(model,6))

        state.keywords['chg-cutoff'] = 10.0
        state.keywords['vdw-cutoff'] = 7.00
        state.keywords['lights'] = ''
        state.keywords['restrainterm'] = ''

        for x in range(0,iter):
            state.dynamics(steps=steps,timestep=1,kw=xtra_kw)
            if not len(state.summary):
                break
            for a in state.summary:
                print(a)
            cmd.load_model(model,'dyna')
            cmd.ending()
            cmd.refresh()
        io.pkl.toFile("realtime.pkl")
        print(" realtime.dyna: terminated after %d steps." % state.counter)

def check(obj='check'):
    global state
    global model

    if not state:
        if not model:
            print(" realtime.reload: please run setup first.")
        else:
            cmd.load_model(model,obj,1)
    else:
        model = state.model
        cmd.load_model(model,obj,1)


def mini(total_step=100,gradient=0.001,interval=100,obj='rt'):

    global state

    if not state:
        print(" realtime.mini: please run setup first.")
    else:
        state.echo = 0

        model = state.model

        print(" realtime.mini: %d atoms total\n" % model.nAtom)

        xtra_kw = []

        xtra_kw.extend(keyword.get_inactive(model,3))
        xtra_kw.extend(keyword.get_restrain_positions(model,2,0,10))
        xtra_kw.extend(keyword.get_restrain_positions(model,5,0.5,1))
        xtra_kw.extend(keyword.get_inactive(model,6))

        state.keywords['chg-cutoff'] = 10.0
        state.keywords['vdw-cutoff'] = 7.00
        state.keywords['lights'] = ''
        state.keywords['restrainterm'] = ''

        iter = total_step/interval
        for x in range(0,iter):
            state.minimize(gradient=gradient,max_iter=interval,kw=xtra_kw)
            cmd.delete(obj)
            cmd.load_model(model,obj,1)
            cmd.refresh()
            if not len(state.summary):
                break
            for a in state.summary:
                print(a)
            if state.summary[-1][7]=='SmallGrad':
                break;
        io.pkl.toFile(model,"realtime.pkl")
        print(" realtime.mini: terminated after %d steps." % state.counter)
