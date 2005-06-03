# commands extensions to PyMOL for batchmin

from pymol import cmd
from chempy.bmin import realtime
import threading

def amin(*arg,**kwarg):
    realtime.assign(arg[0])
    apply(bmin,arg,kwarg)

def bmin(object,iter=500,grad=0.1,interval=100,
            solvation=None):
    realtime.setup(object)
    t = threading.Thread(target=realtime.mini,
                                args=(int(iter),float(grad),int(interval),str(object)),
                                kwargs={
        'solvation':solvation, # turn on GB/SA (cdie)
        'rest_flag':2, # atoms with flag 2 set are restrained
        'fix_flag':3, # atoms with flag 3 set are fixed
        }
                                )
    t.setDaemon(1)
    t.start()

def bmin_sync(object,iter=500,grad=0.1,interval=100,
            solvation=None):
    realtime.setup(object)
    realtime.mini(int(iter),float(grad),int(interval),str(object),
                      solvation=solvation,
                      rest_flag=2,
                      fix_flag=3)
    
def pmin():
    cmd.delete('min')
    realtime.assign('lig')
    cmd.create('min','(lig|prot)')
    bmin('min')
    
cmd.extend('bmin',bmin)
cmd.extend('bmin_sync',bmin_sync)
cmd.extend('amin',amin)
cmd.extend('pmin',pmin)

