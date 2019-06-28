
cmd = __import__("sys").modules["pymol.cmd"]
from pymol import _cmd

def get_progress(reset=0,_self=cmd):
    with _self.lock_api_status:
        return _cmd.get_progress(_self._COb, int(reset))

def ready(_self=cmd): # INTERNAL
    # WARNING: internal routine, subject to change
    return _cmd.ready(_self._COb)
