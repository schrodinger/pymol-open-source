
import threading
import socket
import cPickle
import socket # For gethostbyaddr()
import SocketServer
import sys
import traceback
import copy
import os

# Dictionary Emulator Object Database
#
# After being constructed, DictDBLocal and DictDBClient objects
# work like dictionary objects, but with the following additional methods:
#
# standby() optimizes index file to incorporate changes and purges
#           indexes from memory (indexes will be reread upon next call, if any)
#
# purge() repacks the database and index files to eliminate wasted space
#
# get_info() returns DictDBInfo object containing database statistics
#
# shutdown() terminates database server process (for remote database clients)
#

# PART 1
# local dictionary database class, with a thread-safe API

class DictDBInfo:
   def __init__(self):
      self.used = 0   # bytes containing data
      self.wasted = 0 # bytes wasted (from deleted/overwitten records)

class DictDBLocal:

   __magic__ = r'*8#~'  # 32-bit record stamp for record
   __delete_magic__ = '*8#@'  # 32-bit record stamp for deleted record
   
   def __init__(self,prefix="test_dictdb",bin=1):

      self.index_file = prefix + ".dbi"
      self.data_file = prefix + ".dbf"
      self.bin = bin
      # create lock
      self.lock = threading.RLock()
      # restore indexed into memory
      self._restore()

   def get_info(self):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         result = copy.deepcopy(self.info)
      finally:
         self.lock.release()
      return result

   def has_key(self,key):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         result = self.rec.has_key(key)
      finally:
         self.lock.release()
      return result

   def keys(self):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         result = self.rec.keys()
      finally:
         self.lock.release()
      return result

   def __delitem__(self,key):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         if self.rec.has_key(key):
            del self.rec[key]
            # write Delete magic to data file for recovery
            f=open(self.data_file,'ab')
            f.write(self.__delete_magic__) 
            cPickle.dump(key,f,self.bin) 
            f.close()
            # write blank index to index file
            f=open(self.index_file,'ab')
            cPickle.dump((key,None),f,self.bin)
            f.close()
         else:
            raise KeyError(key)
      finally:
         self.lock.release()
      return result

   def standby(self):
      result = None
      if hasattr(self,'rec'):
         try:
            self.lock.acquire()
            # write index file
            f=open(self.index_file,'wb')
            cPickle.dump(self.info,f,self.bin)
            cPickle.dump(self.rec,f,self.bin)
            f.close()
            del self.rec
            del self.info
         finally:
            self.lock.release()
      return result

   def purge(self):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         tmp_data = self.data_file + '_tmp'
         f = open(self.data_file,'rb')
         g = open(tmp_data,'wb')
         rec = self.rec
         used = 0
         for key in self.rec.keys():
            f.seek(rec[key][0])
            g.write(self.__magic__)
            cPickle.dump(key,g);
            g.write(f.read(rec[key][1]))
            used = used + rec[key][1]
         f.close()
         g.close()
         os.unlink(self.index_file) # delete index file...         
         os.unlink(self.data_file)
         os.rename(tmp_data,self.data_file)
         # write index file
         f=open(self.index_file,'wb')
         cPickle.dump(self.info,f,self.bin)
         cPickle.dump(self.rec,f,self.bin)
         f.close()
         self.info.used = used
         self.info.wasted = 0
      finally:
         self.lock.release()
      self.standby()
      return result

   def shutdown(self): # rebuild index from data file
      try:
         self.lock.acquire()
         sys.exit(0)
      except:
         self.lock.release()
         
   def __getitem__(self,key):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         result = None # default 
         if self.rec.has_key(key):
            f=open(self.data_file)
            f.seek(self.rec[key][0])
            result = cPickle.load(f)
         else:
            raise KeyError(key)
      finally:
         self.lock.release()
      return result

   def _get(self,key): # get with string
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         result = None # default 
         if self.rec.has_key(key):
            f=open(self.data_file)
            f.seek(self.rec[key][0])
            result = f.read(self.rec[key][1])
      finally:
         self.lock.release()
      return result

   def __setitem__(self,key,object):
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         # write data
         f=open(self.data_file,'ab')
         f.write(self.__magic__) # for recovery
         cPickle.dump(key,f,self.bin) # for recovery
         start = f.tell()
         cPickle.dump(object,f,self.bin)
         record_info = (start,f.tell()-start)
         f.close()
         if self.rec.has_key(key):
            self.info.wasted = self.info.wasted + self.rec[key][1]
         self.rec[key] = record_info
         self.info.used = self.info.used + record_info[1]
         # append index file
         f=open(self.index_file,'ab')
         cPickle.dump((key,record_info),f,self.bin)
         f.close()
      finally:
         self.lock.release()
      return result

   def _set(self,key,data_string): # set with string
      result = None
      if not hasattr(self,'rec'):
         self._restore()
      try:
         self.lock.acquire()
         # write data
         f=open(self.data_file,'ab')
         f.write(self.__magic__) # for recovery
         cPickle.dump(key,f,self.bin) # for recovery
         start = f.tell()
         f.write(data_string)
         record_info = (start,f.tell()-start)
         f.close()
         if self.rec.has_key(key):
            self.info.wasted = self.info.wasted + self.rec[key][1]
         self.rec[key] = record_info
         self.info.used = self.info.used + record_info[1]
         # append index file
         f=open(self.index_file,'ab')
         cPickle.dump((key,record_info),f,self.bin)
         f.close()
      finally:
         self.lock.release()
      return result

   def _recover(self): # rebuild index from data file
      result = None
      # need to write recovery routine...
      self.info = DictDBInfo()
      self.rec = {}
      try:
         try:
            self.lock.acquire()
            f=open(self.data_file,'rb')
            # find length of file
            f.seek(0,2)
            eof = f.tell()
            f.seek(0,0)
            if eof:
               print " Warning: Bad index file -- attempting database recovery"
            rec = self.rec
            self_info = self.info
            while f.tell()!=eof:
               chk = f.read(4)
               if chk==self.__magic__: # recover record
                  key = cPickle.load(f)
                  start = f.tell()
                  data = cPickle.load(f)
                  record_info = (start,f.tell()-start)
                  if rec.has_key(key):
                     self_info.wasted = self_info.wasted + rec[key][1]                  
                  rec[key] = record_info
                  self_info.used = self_info.used + record_info[1]
               elif chk==self.__delete_magic__: # delete record
                  key = cPickle.load(f)
                  if rec.has_key(key):
                     self_info.wasted = self_info.wasted + rec[key][1]
                     del rec[key]
               else:
                  raise 'Bad Magic'
         except:
            print " database recovery failed."
            traceback.print_exc()
            sys.exit(1)
         # now write recovered indexes
         f = open(self.index_file,'wb')
         cPickle.dump(self.info,f,self.bin)
         cPickle.dump(self.rec,f,self.bin)
         f.close()
      finally:
         self.lock.release()
      return result
   
   def _restore(self):
      result = None
      try:
         self.lock.acquire()

         # make sure files exist and are writable

         f = open(self.index_file,'ab')
         f.close()

         f = open(self.data_file,'ab')
         f.close()

         # read current indexes

         f = open(self.index_file,'rb')

         # find length of file

         f.seek(0,2)
         eof = f.tell()
         f.seek(0,0)

         # start reading and appending recent changes (if any)

         try:
            self.info = cPickle.load(f)
            self.rec = cPickle.load(f)
            rec = self.rec
            cPickle_load = cPickle.load
            self_info = self.info
            while f.tell()!=eof:
               key, info = cPickle_load(f)
               if info!=None:
                  if rec.has_key(key):
                     self_info.wasted = self_info.wasted + rec[key][1]                     
                  rec[key] = info
                  self_info.used = self_info.used + info[1]
               else:
                  self_info.wasted = self_info.wasted + rec[key][1]
                  del rec[key]
            f.close()
         except EOFError:
            # if error reading indexes, then do datafile recovery
            f.close()
            self._recover()
      finally:
         self.lock.release()
      return result

# PART 2
# Database client

class DictDBClient:

   def __init__(self,host='localhost',port=8000):
      self.host=host
      self.port=port

   def _remote_call(self,meth,args,kwds):
      result = None
      s=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
      s.connect((self.host,self.port))
      send = s.makefile('w')
      recv = s.makefile('r')
      cPickle.dump(meth,send,1) # binary by default
      cPickle.dump(args,send,1)
      cPickle.dump(kwds,send,1)
      send.close()
      result = cPickle.load(recv)
      recv.close()
      s.close()      
      return result

   def __getitem__(self,key):
      if self._remote_call('has_key',(key,),{}):
         return cPickle.loads(self._remote_call('_get',(key,),{}))
      else:
         raise KeyError(key)

   def __setitem__(self,key,object):
      return self._remote_call('_set',(key,cPickle.dumps(object)),{})

   def __delitem__(self,key):
      if self._remote_call('has_key',(key,),{}):
         return self._remote_call('__delitem__',(key,),{})
      else:
         raise KeyError(key)

   def standby(self,key):
      return self._remote_call('standby',(),{})

   def purge(self):
      return self._remote_call('purge',(),{})

   def keys(self):
      return self._remote_call('keys',(),{})

   def has_key(self,key):
      return self._remote_call('has_key',(key,),{})

   def get_info(self):
      return self._remote_call('get_info',(),{})

def server_test(port = 8000):
   print 'Testing DictDBServer on port',str(port)
   fp = DictDBServer(port=port) # socket servers don't terminate

# PART 3
# Database Socket Server

class DictDBServer:
   def __init__(self,prefix='',port=''):
      
      server_address = ('', port)

      ddbs = _DictDBServer(server_address, DictDBRequestHandler)  

      # assign dict database to this server
      ddbs.dictdb = DictDBLocal(prefix)

      # now serve requests forever
      ddbs.serve_forever()
      
class _DictDBServer(SocketServer.ThreadingTCPServer):

    def server_bind(self):
        """Override server_bind to store the server name."""
        SocketServer.ThreadingTCPServer.server_bind(self)
        host, port = self.socket.getsockname()
        if not host or host == '0.0.0.0':
            host = socket.gethostname()
        hostname, hostnames, hostaddrs = socket.gethostbyaddr(host)
        if '.' not in hostname:
            for host in hostnames:
                if '.' in host:
                    hostname = host
                    break
        self.server_name = hostname
        self.server_port = port

class DictDBRequestHandler(SocketServer.StreamRequestHandler):

    def handle(self):
        # get method name from client
        method = cPickle.load(self.rfile)
        # get arguments from client
        args = cPickle.load(self.rfile)
        kw = cPickle.load(self.rfile)
        # get method pointer
        meth_obj = getattr(self.server.dictdb,method)
        print method,args,kw
        # call method and return result
        cPickle.dump(apply(meth_obj,args,kw),self.wfile,1) # binary by default

def client_test(host,port=8000):

   print 'Testing Client with Server on port',str(port)

   tdb = DictDBClient(port=port)
   tdb['test']='hello'
   print tdb['test']
   try:
      print tdb['nonexistent']
   except KeyError:
      print " key error 1 as expected"
   try:
      del tdb['nonexistent']
   except KeyError:
      print " key error 2 as expected"
   tdb['extra'] = 'hi'
   print tdb.has_key('extra')
   print tdb.keys()

if __name__=='__main__':
   import os
   
   print '***Testing DictDB***:'

   ddb = DictDBLocal()
   print ddb.keys()
   ddb['test']='some data object'
   print ddb['test']
   ddb['another']='another data object'
   print ddb['another']
   ddb['test']='some updated data object'   
   print ddb['test']
   del ddb

   ddb = DictDBLocal()
   ddb['whoa']='whoa data object'
   print ddb['test'] 
   ddb['dude']='dude data object'
   print ddb['dude']
   ddb['number'] = 9999
   print ddb['number']
   ddb.standby()
   del ddb

   ddb = DictDBLocal()
   print 'Current keys:'+str(ddb.keys())

   print ddb['test']
   print ddb['another']

   print '***Testing raw string methods:'
   print cPickle.loads(ddb._get('test'))
   print cPickle.loads(ddb._get('another'))

   ddb._set('some_key',cPickle.dumps("raw test1"))
   print ddb['some_key']

   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   del ddb
   
   os.unlink("test_dictdb.dbi")
   ddb = DictDBLocal()
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   del ddb

   os.unlink("test_dictdb.dbi")
   ddb = DictDBLocal()
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   ddb.standby()
   del ddb

   ddb = DictDBLocal()
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   del ddb

   ddb = DictDBLocal()
   ddb.standby()
   ddb['reactivated'] = 1234
   print ddb.has_key('reactivated')
   ddb.standby()
   ddb['reactivated'] = 495
   print ddb['reactivated']
   try:
      print ddb['nonexistent']
   except KeyError:
      print " key error 1 as expected"
   try:
      del ddb['nonexistent']
   except KeyError:
      print " key error 2 as expected"
   
   print "purge test:"
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   ddb.purge()
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   del ddb

   ddb = DictDBLocal()
   info = ddb.get_info()
   print info.used,info.wasted,info.used-info.wasted
   del ddb


   
   

