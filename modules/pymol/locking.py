
import thread
import threading

import pymol
import cmd

from cmd import fb_module, fb_mask, fb_action

import _cmd

# the following lock is used by both C and Python to insure that no more than
# one active thread enters PyMOL at a given time. 

lock_api = pymol.lock_api
lock_api_c = pymol.lock_api_c
lock_api_status = pymol.lock_api_status
lock_api_glut = pymol.lock_api_glut


# WARNING: internal routines, subject to change      
def lock_c(_self=cmd): 
    lock_api_c.acquire(1)

def unlock_c(_self=cmd):
    lock_api_c.release()

def lock_status_attempt(_self=cmd):
    return lock_api_status.acquire(0)

def lock_status(_self=cmd): 
    lock_api_status.acquire(1)

def unlock_status(_self=cmd):
    lock_api_status.release()

def lock_glut(_self=cmd): 
    lock_api_glut.acquire(1)

def unlock_glut(_self=cmd):
    lock_api_glut.release()

def lock_without_glut(_self=cmd):
    try:
        lock_glut()
        lock(_self)
    finally:
        unlock_glut()

def lock(_self=cmd): # INTERNAL -- API lock
#      print " lock: acquiring as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
    if not lock_api.acquire(0):
        w = 0.001
        while 1:
#            print " lock: ... as 0x%x"%thread.get_ident(),(thread.get_ident() == pymol.glutThread)
            e = threading.Event() 
            e.wait(w)  
            del e
            if lock_api.acquire(0):
                break
            if w<0.1:
                w = w * 2 # wait twice as long each time until flushed
#      print "lock: acquired by 0x%x"%thread.get_ident()

def lock_attempt(_self=cmd): # INTERNAL
    return lock_api.acquire(blocking=0)

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
        lock_api.release()
    #         print "lock: released by 0x%x (glut)"%thread.get_ident()
        if result==None: # don't flush if we have an incipient error (negative input)
            _cmd.flush_now(_self._COb)
        elif cmd.is_ok(result):
            _cmd.flush_now(_self._COb)
    else:
    #         print "lock: released by 0x%x (not glut), waiting queue"%thread.get_ident()
        lock_api.release()
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
                    if _self._feedback(fb_module.cmd,fb_mask.debugging):
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
