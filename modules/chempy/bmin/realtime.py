# pymol

from pymol import cmd
from chempy import io
from chempy import feedback
from chempy.bmin.state import State

import threading
import traceback
import os

state = None
model = None

def assign(sele,preserve=0):

    from molobj import MolObj
    from typer import Typer,Rules
    import rules

    result = 1
    
    state = State()

    model = cmd.get_model(sele)
    # now assign atom types

    ruleSet = Rules()

    ruleSet.fromList(rules.mmff_types)
    ruleSet.mappingFromList(rules.mmff_mapping)
    
    mobj = MolObj()
    mobj.fromChemPyModel(model)

    typed = Typer(molObj = mobj)
    
    print " realtime: assigning atom types"
    typed.applyRules(ruleSet)

    c = 0
    for a in typed.getTypes():
        at = model.atom[c]
        if (at.text_type == '??') or (not preserve):
            if a==-99:
                print " warning: unable to assign atom type to atom %d"%c
                result = 0
            else:
                cmd.alter("((%s) and (index %s))" % (sele,at.index),
                             "numeric_type ='%s'" % a)
                if feedback['tinker']:
                    print " "+str(__name__)+': %s is a %s' % (at.name,a)
                at.numeric_type = a
        c = c + 1

    sm = 0
    for a in model.atom:
        a.resi = str(a.resi_number)
        sm = sm + a.partial_charge

    return result

def setup(sele,preserve=0):
    
    global state
    global model
    
    state = State()
    model = cmd.get_model(sele)

    sm = 0
    for a in model.atom:
        a.resi = str(a.resi_number)
        sm = sm + a.partial_charge

    state.load_model(model)
    return 1

def check(obj='check'):
    global state
    global model
    
    if not state:
        if not model:
            print " realtime.reload: please run setup first."
        else:
            cmd.load_model(model,obj,1)
    else:
        model = state.model
        cmd.load_model(model,obj,1)

def mini(total_steps=500,
            gradient=0.001,
            interval=100,
            object='rt',
            fix_flag=None,
            rest_flag=None,
            solvation=None,
            finish=None):

    global state
    if not state:
        print " realtime.mini: please run setup first..."
    else:
        model = state.model
        print " realtime.mini: %d atoms total\n" % model.nAtom
        try:
            while total_steps>0:
                total_steps = total_steps - interval
                state.minimize(max_iter=interval,
                                    fix_flag=fix_flag,
                                    rest_flag=rest_flag,
                                    solvation=solvation)
                cmd.delete(object)
                cmd.load_model(state.model,object,1)
                cmd.refresh()
                if finish!=None:
                    apply(finish[0],finish[1],finish[2])
        except:
            cmd.load_model(state.model,'ref')
            traceback.print_exc()
        print " realtime.mini: complete."

def mini_threaded(*args,**kwargs):
    t = threading.Thread(target=mini,
                                args=args,
                                kwargs=kwargs)
    t.setDaemon(1)
    t.start()
    



