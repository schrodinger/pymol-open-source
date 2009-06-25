
import thread
import threading

import pymol
import cmd

from cmd import fb_module, fb_mask, fb_action, fb_debug

import _cmd

# WARNING: internal routines, subject to change      
def lock_c(_self=cmd): 
    _self.lock_api_c.acquire(1)

def unlock_c(_self=cmd):
    _self.lock_api_c.release()

def lock_data(_self=cmd):
    _self.lock_api_data.acquire(1)

def unlock_data(_self=cmd):
    _self.lock_api_data.release()
    
def lock_status_attempt(_self=cmd):
    return _self.lock_api_status.acquire(0)

def lock_status(_self=cmd): 
    _self.lock_api_status.acquire(1)

def unlock_status(_self=cmd):
    _self.lock_api_status.release()

def lock_glut(_self=cmd): 
    _self.lock_api_glut.acquire(1)

def unlock_glut(_self=cmd):
    _self.lock_api_glut.release()

def lock_without_glut(_self=cmd):
    try:
        _self.lock_glut(_self)
        _self.lock(_self)
    finally:
        _self.unlock_glut(_self)

def lock(_self=cmd): # INTERNAL -- API lock
#      print " lock: acquiring as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
    if not _self.lock_api.acquire(0):
        w = 0.0001
        while 1:
#            print " lock: ... as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
            e = threading.Event() 
            e.wait(w)  
            del e
            if _self.lock_api.acquire(0):
                break
            if w<0.05:
                w = w * 2 # wait twice as long each time until flushed
            else: # we're not getting lucky, so block for real
                _self.lock_api.acquire(1)
                break
#      print "lock: acquired by 0x%x"%thread.get_ident()

def lock_attempt(_self=cmd): # INTERNAL
    return _self.lock_api.acquire(blocking=0)

def block_flush(_self=cmd):
    lock(_self)
    _self.lock_api_allow_flush = 0
    unlock(None,_self)

def unblock_flush(_self=cmd):
    lock(_self)
    _self.lock_api_allow_flush = 1
    unlock(None,_self)
    
def unlock(result=None,_self=cmd): # INTERNAL
    if (thread.get_ident() == pymol.glutThread):
        if _self.reaper:
            try:
                if not _self.reaper.isAlive():
                    if pymol.invocation.options.no_gui:
                        _cmd.quit(_self._COb) # TO FIX
                    else:
                        _self.reaper = None
            except:
                pass
        _self.lock_api.release()
    #         print "lock: released by 0x%x (glut)"%thread.get_ident()
        if _self.lock_api_allow_flush:
            if result==None: # don't flush if we have an incipient error (negative input)
                _cmd.flush_now(_self._COb)
            elif _self.is_ok(result):
                
                _cmd.flush_now(_self._COb)
    else:
    #         print "lock: released by 0x%x (not glut), waiting queue"%thread.get_ident()
        _self.lock_api.release()
        if _cmd.wait_queue(_self._COb): # commands waiting to be executed?
            e = threading.Event() # abdicate control for a 100 usec for quick tasks
            e.wait(0.0001)
            del e
            # then give PyMOL increasingly longer intervals to get its work done...
            w = 0.0005  # NOTE: affects API perf. for "do" and delayed-exec
            while _cmd.wait_queue(_self._COb):
                e = threading.Event() # abdicate control for a 100 usec for quick tasks
                e.wait(w)
                del e
                if w > 0.1: # wait up 0.2 sec max for PyMOL to flush queue
                    if _self._feedback(fb_module.cmd,fb_mask.debugging,_self):
                        fb_debug.write("Debug: avoiding possible dead-lock?\n")
   #                      print "dead locked as 0x%x"%thread.get_ident()
                    break
                w = w * 2 # wait twice as long each time until flushed
            
def is_glut_thread(): # internal
    if thread.get_ident() == pymol.glutThread:
        return 1
    else:
        return 0

def setup_global_locks(): # INTERNAL, OBSOLETE?
    # WARNING: internal routine, subject to change
    pass

def interrupt(_self=cmd): # asynch -- no locking!
    _cmd.interrupt(_self._COb,1)
    return None
