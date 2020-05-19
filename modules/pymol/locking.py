
import sys
import _thread as thread
import threading

cmd = sys.modules["pymol.cmd"]

from pymol import _cmd

# WARNING: internal routines, subject to change
def lock_without_glut(_self=cmd):
    with _self.lock_api_glut:
        _self.lock(_self)

class LockCM(object):
    '''
    API lock context manager
    '''
    def __init__(self, _self=cmd):
        self.cmd = _self._weakrefproxy
    def __enter__(self):
        lock(self.cmd)
    def __exit__(self, type, value, traceback):
        unlock(None if type is None else -1, self.cmd)

def lock(_self=cmd): # INTERNAL -- API lock
    return _self.lock_api.acquire()

def lock_attempt(_self=cmd): # INTERNAL
    return _self.lock_api.acquire(blocking=0)

def block_flush(_self=cmd):
    with _self.lockcm:
        _self.lock_api_allow_flush = 0

def unblock_flush(_self=cmd):
    with _self.lockcm:
        _self.lock_api_allow_flush = 1

def unlock(result=None,_self=cmd): # INTERNAL
    '''
    Release the API lock and flush the command queue
    '''
    thread_owns_gui = _self.is_gui_thread()

    if thread_owns_gui and _self.reaper is not None:
        # TODO what's the use case? How can the finish_launching() thread
        # die and what's the expected behavior? (see PYMOL-3247)
            try:
                if not _self.reaper.is_alive():
                    if _self._pymol.invocation.options.no_gui:
                        _cmd.quit(_self._COb) # TO FIX
                    else:
                        _self.reaper = None
            except:
                pass

    _self.lock_api.release()

    # don't flush if we have an incipient error (negative input)
    if _self.is_error(result):
        return

    if thread_owns_gui:
        if _self.lock_api_allow_flush:
            _cmd.flush_now(_self._COb)
    else:
        # TODO: what's the logic here? We assume we have to wait for
        # cmd.do() calls being executed by GUI thread, but only for 200ms?
        # What's the concrete use case and would there be a better (more
        # predictable) logic for it? Isn't this an invitation for race
        # conditions, e.g. other threads could call cmd.do() and make us
        # wait for it? (see PYMOL-3248)
        # NOTE: affects API perf. for "do" and delayed-exec
        w = 0.0005
        while w < 0.1 and _cmd.wait_queue(_self._COb):
            threading.Event().wait(w)
            w *= 2

def is_gui_thread(_self=cmd): # internal
    '''
    Return true if the current thread is the GUI thread (e.g. GLUT or Qt)
    or if there is no GUI thread. False otherwise.
    '''
    gui_ident = _self._pymol.glutThread
    return gui_ident is None or gui_ident == thread.get_ident()

def interrupt(_self=cmd): # asynch -- no locking!
    _cmd.interrupt(_self._COb,1)
