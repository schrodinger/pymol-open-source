# simple proof-of-concept exercise in PyMOL synchronization

# NOT FOR PRODUCTION USE! HIGHLY INSECURE!

# You need to launch two PyMOL simultaneously as follows...

# MASTER/SLAVE (on same machine):

# pymol syncmol.py -- recv 8000
# pymol syncmol.py -- send 8000

# MUTUAL SYNCHRONIZATION (on same machine):

# pymol syncmol.py -- recv 8001 send 8002
# pymol syncmol.py -- send 8001 recv 8002

# MUTUAL SYNCHRONIZATION (over a lan):

# on host1: pymol syncmol.py -- recv 8000 send host2:8000
# on host2: pymol syncmol.py -- recv 8000 send host1:8000

import threading
import socket
import cPickle
import socket # For gethostbyaddr()
import SocketServer
import sys
import traceback
import copy
import os
import time
import Queue

    
class PyMOLWriter: # this class transmits

    def __init__(self, pymol, host='localhost', port=8000):
        host = str(host)
        port = int(port)
        self.host = host
        self.port = port
        self.sock = None
        self.cmd = pymol.cmd
        self.fifo = Queue.Queue(0)
        cmd = self.cmd

        print " syncmol: writing to %s:%d"%(host,port)
        pymol.cmd.log_open(self.fifo)
        
        last_view = None
        last_frame = 0
        while 1:
            time.sleep(0.1) # update 10x a second
            view = cmd.get_view(output=4)
            try:
                cmd.lock()
                if hasattr(cmd,"_last_view"):
                    last_view = cmd._last_view
                    del cmd._last_view
            finally:
                cmd.unlock(-1)
            deferred = 0
            if not self.fifo.empty():
                if view != last_view:
                    self._remote_call("do",("_ cmd.set('defer_updates')",),{'log':0})
                    deferred = 1
                do_list = []
                while not self.fifo.empty():
                    do_list.append(self.fifo.get())
                self._remote_call("do",(do_list,),{'log':0})
            if view != last_view:
                self._remote_call("remote_set_view",(view,))
                last_view = view
                if deferred:
                    self._remote_call("do",("_ cmd.unset('defer_updates')",),{'log':0})
            if not cmd.get_movie_playing():
                frame = int(cmd.get("frame"))
                if last_frame != frame:
                    self._remote_call("frame",(frame,))
                    last_frame = frame
            else:
                last_frame = None
            
                
    def _remote_call(self,meth,args=(),kwds={}):
        result = None
        while self.sock == None:
            try:
                self.sock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
                self.sock.connect((self.host,self.port))
                self.send = self.sock.makefile('w')
                self.recv = self.sock.makefile('r')
            except:
                self.sock=None
                print "retrying..."
                time.sleep(1)
        cPickle.dump(meth,self.send,1) # binary by default
        cPickle.dump(args,self.send,1)
        cPickle.dump(kwds,self.send,1)
        self.send.flush()
        result = cPickle.load(self.recv)
        return result

def remote_set_view(view,_self=cmd):
    try:
        _self.lock()
        _self._last_view = view
        _self.set_view(view)
    finally:
        _self.unlock(-1)
    
class PyMOLReader: # this class receives

    def __init__(self,pymol,port):

        sys.setcheckinterval(0)
        
        server_address = ('', port)

        ddbs = _PyMOLReader(server_address, _PyMOLRequestHandler)  

        # bind pymol instance to the reader

        pymol.cmd.remote_set_view = remote_set_view
        ddbs.cmd = pymol.cmd

        print " syncmol: reading from port %d"%(port)
        # now serve requests forever
        ddbs.keep_alive = 1
        while ddbs.keep_alive:
            ddbs.handle_request()
        
class _PyMOLReader(SocketServer.ThreadingTCPServer):

     def server_bind(self):
          """Override server_bind to store the server name."""
          SocketServer.ThreadingTCPServer.server_bind(self)
          host, port = self.socket.getsockname()
          if not host or host == '0.0.0.0':
                host = socket.gethostname()
          try:
              hostname, hostnames, hostaddrs = socket.gethostbyaddr(host)
              if '.' not in hostname:
                    for host in hostnames:
                         if '.' in host:
                              hostname = host
                              break
          except:
              hostname = 'localhost'
          self.server_name = hostname
          self.server_port = port

class _PyMOLRequestHandler(SocketServer.StreamRequestHandler):

     def handle(self):
         while self.server.keep_alive:
             # get method name from client

             try:
                 method = cPickle.load(self.rfile)
             except EOFError,socket.error:
                 break

             if method == 'shutdown':
                 self.server.keep_alive = 0
                 
             # get arguments from client
             args = cPickle.load(self.rfile)
             kw = cPickle.load(self.rfile)

             # get cmd method pointer
             meth_obj = getattr(self.server.cmd,method)

             # call method and return result
             cPickle.dump(apply(meth_obj,args,kw),self.wfile,1) # binary by default
             self.wfile.flush()
             
if __name__=='pymol':
    import os
    import pymol
    sys.argv.reverse()
    sys.argv.pop()
    while len(sys.argv):
        tok = sys.argv.pop()
        if tok == 'recv':
            port = int(sys.argv.pop())
            _stdin_reader_thread = threading.Thread(target=PyMOLReader,
                                                    args=(pymol,port))
            _stdin_reader_thread.setDaemon(1)
            _stdin_reader_thread.start()
        elif tok == 'send':
            addr = string.split(sys.argv.pop(),':')
            if len(addr)==1:
                host = 'localhost'
                port = int(addr[0])
            else:
                host = addr[0]
                port = addr[1]
            _stdin_reader_thread = threading.Thread(target=PyMOLWriter,
                                                    args=(pymol,host,port))
            _stdin_reader_thread.setDaemon(1)
            _stdin_reader_thread.start()

            
    
