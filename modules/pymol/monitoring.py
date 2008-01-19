
import cmd
import _cmd
import sys

def get_progress(reset=0,_self=cmd):
    r = -1.0
    try:
        _self.lock_status(_self)
        r = _cmd.get_progress(_self._COb,int(reset))
    finally:
        _self.unlock_status(_self)
    return r

def check_redundant_open(file,_self=cmd):
    found = 0
    for a in _self._pymol.invocation.options.deferred:
        if a == file:
            found = 1
            break
    for a in _self._pymol.invocation._argv: 
        if a == file:
            found = 1
            break;
    return found

def ready(_self=cmd): # INTERNAL
    # WARNING: internal routine, subject to change
    return _cmd.ready(_self._COb)
