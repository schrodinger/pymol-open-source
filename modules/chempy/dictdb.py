#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------

import threading
import socket
import socket # For gethostbyaddr()
import sys
if True:
    import pickle as cPickle
    import socketserver as SocketServer
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
# shutdown() terminates the database server (only useful for remote clients)
#
# reset() clears database

# PART 1
# Database Engine, with a thread-safe API

class DictDBInfo:
    def __init__(self):
        self.used = 0   # bytes containing data
        self.wasted = 0 # bytes wasted (from deleted/overwitten records)

class DictDBLocal:

    __magic__ = b'*8#~'  # 32-bit record stamp for record
    __delete_magic__ = b'*8#@'  # 32-bit record stamp for deleted record


    def __init__(self,prefix,bin=1,read_only=0):

        # store important information
        self.index_file = prefix + ".dbi"  # database information (can be reconstructed from .dbf)
        self.data_file = prefix + ".dbf"   # database data
        self.bin = bin
        self.changed = 0
        self.read_only = read_only

        # create lock
        self.lock = threading.RLock() # NOTE: recursive for convenience

        # restore indexed into memory
        self._restore()


    def get_info(self):
        result = None

        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()

            # generate independent copy for calling thread
            result = copy.deepcopy(self.info)
        finally:
            self.lock.release()
        return result


    def has_key(self,key):
        result = None

        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            result = key in self.rec
        finally:
            self.lock.release()
        return result

    __contains__ = has_key

    def keys(self):
        result = None

        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()

            # get list of keys (thread safe...result is a new object)
            result = list(self.rec.keys())
        finally:
            self.lock.release()
        return result


    def __delitem__(self,key):
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            if key in self.rec:

                # account for space
                if key in self.rec:
                    self.info.wasted = self.info.wasted + self.rec[key][1]

                # delete record
                del self.rec[key]

                if not self.read_only:
                    # write delete magic to data file for recovery
                    f=open(self.data_file,'ab')
                    f.write(self.__delete_magic__)
                    cPickle.dump(key,f,self.bin)
                    f.close()

                    # write blank index to index file
                    f=open(self.index_file,'ab')
                    cPickle.dump((key,None),f,self.bin)
                    f.close()

                # note change
                self.changed = 1
            else:
                raise KeyError(key)
        finally:
            self.lock.release()
        return result


    def standby(self):
        result = None

        # only do something if dictionary is already in RAM
        if hasattr(self,'rec') and not self.read_only:
            try:
                self.lock.acquire()

                # write index file
                f=open(self.index_file,'wb')
                cPickle.dump(self.info,f,self.bin)
                cPickle.dump(self.rec,f,self.bin)
                f.close()

                # free memory
                del self.rec
                del self.info
            finally:
                self.lock.release()
        return result

    def reset(self):
        result = None
        try:
            self.lock.acquire()

            if not self.read_only:
            # make sure files exist and are writable

                f = open(self.index_file,'wb')
                f.close()

                f = open(self.data_file,'wb')
                f.close()

            # new database information

            self.info = DictDBInfo()
            self.rec = {}
            self.changed = 1
            self.standby() # write out blank indexes (important)

        finally:
            self.lock.release()
        return result

    def purge(self):
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            if not self.read_only:
                # open source and temporary data files
                tmp_data = self.data_file + '_tmp'
                f = open(self.data_file,'rb')
                g = open(tmp_data,'wb')
                rec = self.rec
                used = 0

                # iterate through records, copying only those which are extant
                for key in list(self.rec.keys()):
                    f.seek(rec[key][0])
                    g.write(self.__magic__)
                    cPickle.dump(key,g,self.bin);
                    start = g.tell()
                    g.write(f.read(rec[key][1]))
                    used = used + rec[key][1]
                    rec[key]=(start,rec[key][1]) # generate new record entry
                f.close()
                g.close()

                # now perform the switch-over
                os.unlink(self.index_file) # delete old index file...
                os.unlink(self.data_file)  # delete old data file
                os.rename(tmp_data,self.data_file) # move new data file over old

                # update database information
                self.info.used = used
                self.info.wasted = 0

                # write new index and information file
                f=open(self.index_file,'wb')
                cPickle.dump(self.info,f,self.bin)
                cPickle.dump(self.rec,f,self.bin)
                f.close()
        finally:
            self.lock.release()
        return result

    def shutdown(self): # dummy
        return None

    def __getitem__(self,key): # get with object
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            if key in self.rec:

                # locate and retrieve object
                f=open(self.data_file,'rb')
                f.seek(self.rec[key][0])
                result = cPickle.load(f)
                f.close()
            else:
                raise KeyError(key)
        finally:
            self.lock.release()
        return result


    def _get(self,key): # get with string (for remote connections)
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            if key in self.rec:

                # locate and retrieve string
                f=open(self.data_file,'rb')
                f.seek(self.rec[key][0])
                result = f.read(self.rec[key][1])
                f.close()
        finally:
            self.lock.release()
        return result


    def __setitem__(self,key,object):
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()

            # append data onto data file
            f=open(self.data_file,'ab')
            f.write(self.__magic__) # for recovery
            cPickle.dump(key,f,self.bin) # for recovery
            start = f.tell()
            cPickle.dump(object,f,self.bin)
            record_info = (start,f.tell()-start)
            f.close()

            # account for space (if replacing)
            if key in self.rec:
                self.info.wasted = self.info.wasted + self.rec[key][1]

            # update record info
            self.rec[key] = record_info

            # account for space
            self.info.used = self.info.used + record_info[1]

            # append new record  onto index file
            f=open(self.index_file,'ab')
            cPickle.dump((key,record_info),f,self.bin)
            f.close()

            # note change
            self.changed = 1
        finally:
            self.lock.release()
        return result


    def _set(self,key,data_string): # set with string (for remote connections)
        result = None
        # restore dictionary (if nec.)
        if not hasattr(self,'rec'):
            self._restore()
        try:
            self.lock.acquire()
            self.changed = 1

            if not self.read_only:
                # append data onto data file
                f=open(self.data_file,'ab')
                f.write(self.__magic__) # for recovery
                cPickle.dump(key,f,self.bin) # for recovery
                start = f.tell()
                f.write(data_string)
                record_info = (start,f.tell()-start)
                f.close()

            # account for space (if replacing)
            if key in self.rec:
                self.info.wasted = self.info.wasted + self.rec[key][1]

            # update record info
            self.rec[key] = record_info

            # account for space
            self.info.used = self.info.used + record_info[1]

            if not self.read_only:
                # append new record onto index file
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

                # create locals for better performance
                rec = self.rec
                self_info = self.info
                while f.tell()!=eof:
                    chk = f.read(4)
                    if chk==self.__magic__: # recover extant record
                        key = cPickle.load(f)
                        start = f.tell()
                        data = cPickle.load(f)
                        record_info = (start,f.tell()-start)

                        # account for space (if replacing)
                        if key in rec:
                            self_info.wasted = self_info.wasted + rec[key][1]
                        rec[key] = record_info

                        # account for psace
                        self_info.used = self_info.used + record_info[1]
                    elif chk==self.__delete_magic__: # delete already recovered record
                        key = cPickle.load(f)
                        if key in rec:
                            # account for space
                            self_info.wasted = self_info.wasted + rec[key][1]
                            del rec[key]
                    else:
                        raise RuntimeError('Bad Magic')
            except:
                print(" dictdb error: database recovery failed.")
                traceback.print_exc()
                sys.exit(1)
            if not self.read_only:
                # write recovered indexes to a new index file
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

                # use locals for better performance
                rec = self.rec
                cPickle_load = cPickle.load
                self_info = self.info

                while f.tell()!=eof:
                    key, info = cPickle_load(f)
                    if info is not None:

                        # account for space (if replacing)
                        if key in rec:
                            self_info.wasted = self_info.wasted + rec[key][1]

                        # set record infor
                        rec[key] = info

                        # account for space
                        self_info.used = self_info.used + info[1]
                    else:
                        # account for space (if deleting object)
                        self_info.wasted = self_info.wasted + rec[key][1]
                        del rec[key]
                f.close()
            except EOFError:
                # if error occurs when reading indexes, then do a datafile-based recovery
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
        self.sock = None

    def _remote_call(self,meth,args,kwds):
        result = None
        if self.sock is None:
            self.sock=socket.socket(socket.AF_INET,socket.SOCK_STREAM)
            self.sock.connect((self.host,self.port))
            self.send = self.sock.makefile('w')
            self.recv = self.sock.makefile('r')
        cPickle.dump(meth,self.send,1) # binary by default
        cPickle.dump(args,self.send,1)
        cPickle.dump(kwds,self.send,1)
        self.send.flush()
        result = cPickle.load(self.recv)
        return result

    def __getitem__(self,key):
        if self._remote_call('has_key',(key,),{}):
            return cPickle.loads(self._remote_call('_get',(key,),{}))
        else:
            raise KeyError(key)

    def __setitem__(self,key,object):
        self.changed = 1
        return self._remote_call('_set',(key,cPickle.dumps(object)),{})


    def __delitem__(self,key):
        if self._remote_call('has_key',(key,),{}):
            return self._remote_call('__delitem__',(key,),{})
        else:
            raise KeyError(key)

    def standby(self):
        return self._remote_call('standby',(),{})

    def shutdown(self):
        try:
            # multiple connections are sometimes required...
            self._remote_call('shutdown',(),{})
            self.sock=None
            self._remote_call('shutdown',(),{})
            self.sock=None
            self._remote_call('shutdown',(),{})
            self.sock=None
            self._remote_call('shutdown',(),{})
        except:
            pass
        return None

    def purge(self):
        return self._remote_call('purge',(),{})

    def reset(self):
        return self._remote_call('reset',(),{})

    def keys(self):
        return self._remote_call('keys',(),{})

    def has_key(self,key):
        return self._remote_call('has_key',(key,),{})

    __contains__ = has_key

    def get_info(self):
        return self._remote_call('get_info',(),{})

# PART 3
# Database Socket Server

class DictDBServer:
    def __init__(self,prefix,port=''):

        sys.setcheckinterval(0)

        server_address = ('', port)

        ddbs = _DictDBServer(server_address, DictDBRequestHandler)

        # assign dict database to this server
        ddbs.dictdb = DictDBLocal(prefix)

        # now serve requests forever
        ddbs.keep_alive = 1
        while ddbs.keep_alive:
            ddbs.handle_request()

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
         while self.server.keep_alive:
             # get method name from client

             try:
                 method = cPickle.load(self.rfile)
             except (EOFError, socket.error):
                 break

             if method == 'shutdown':
                 self.server.keep_alive = 0

             # get arguments from client
             args = cPickle.load(self.rfile)
             kw = cPickle.load(self.rfile)

             # get method pointer
             meth_obj = getattr(self.server.dictdb,method)
#          print method,args,kw

             # call method and return result
             cPickle.dump(meth_obj(*args, **kw),self.wfile,1) # binary by default
             self.wfile.flush()

def server_test(port = 8000,prefix='test_dictdb'):
    print('Testing DictDBServer on port',str(port))
    fp = DictDBServer(prefix,port=port) # socket servers don't terminate

def client_test(host,port=8000):

    print('Testing Client with Server on port',str(port))

    tdb = DictDBClient(port=port)
    tdb['test']='hello'
    print(tdb['test'])
    try:
        print(tdb['nonexistent'])
    except KeyError:
        print(" key error 1 as expected")
    try:
        del tdb['nonexistent']
    except KeyError:
        print(" key error 2 as expected")
    tdb['extra'] = 'hi'
    print('extra' in tdb)
    print(list(tdb.keys()))

if __name__=='__main__':
    import os

    print('***Testing DictDB***:')

    ddb = DictDBLocal('test_dictdb')
    print(list(ddb.keys()))
    ddb['test']='some data object'
    print(ddb['test'])
    ddb['another']='another data object'
    print(ddb['another'])
    ddb['test']='some updated data object'
    print(ddb['test'])
    del ddb

    ddb = DictDBLocal('test_dictdb')
    ddb['whoa']='whoa data object'
    print(ddb['test'])
    ddb['dude']='dude data object'
    print(ddb['dude'])
    ddb['number'] = 9999
    print(ddb['number'])
    ddb.standby()
    del ddb

    ddb = DictDBLocal('test_dictdb')
    print('Current keys:'+str(list(ddb.keys())))

    print(ddb['test'])
    print(ddb['another'])

    print('***Testing raw string methods:')
    print(cPickle.loads(ddb._get('test')))
    print(cPickle.loads(ddb._get('another')))

    ddb._set('some_key',cPickle.dumps("raw test1"))
    print(ddb['some_key'])

    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    del ddb

    os.unlink("test_dictdb.dbi")
    ddb = DictDBLocal('test_dictdb')
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    del ddb

    os.unlink("test_dictdb.dbi")
    ddb = DictDBLocal('test_dictdb')
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    ddb.standby()
    del ddb

    ddb = DictDBLocal('test_dictdb')
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    del ddb

    ddb = DictDBLocal('test_dictdb')
    ddb.standby()
    ddb['reactivated'] = 1234
    print('reactivated' in ddb)
    ddb.standby()
    ddb['reactivated'] = 495
    print(ddb['reactivated'])
    try:
        print(ddb['nonexistent'])
    except KeyError:
        print(" key error 1 as expected")
    try:
        del ddb['nonexistent']
    except KeyError:
        print(" key error 2 as expected")

    print("purge test:")
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    ddb.purge()
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    del ddb

    ddb = DictDBLocal('test_dictdb')
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)

    ddb['hello']=10
    ddb[123] = 'bacon'
    ddb['green']= 'color'

    print(ddb['hello'])

    print(list(ddb.keys()))

    print('green' in ddb)
    print('blue' in ddb)

    del ddb['hello']

    try:
        print(ddb['hello'])
    except KeyError:
        print(' got expected key error 1')

    try:
        del ddb['hello']
    except KeyError:
        print(' got expected key error 2')

    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    ddb.purge()
    info = ddb.get_info()
    print(info.used,info.wasted,info.used-info.wasted)
    print(ddb['green'])

    import random

    ddc=DictDBLocal('test_dictdb')

    print("loading...")

    lst = []
    for a in range(1,250):
        ddc[a]=str(a)
        lst.append([random.random(),a])
        ddc[str(a)] = a

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print('verifying...')
    for a in range(1,250):
        if ddc[a]!=str(a):
            raise RuntimeError

    for a in range(1,250):
        del ddc[str(a)]

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print('verifying...')
    for a in range(1,250):
        if ddc[a]!=str(a):
            raise RuntimeError

    del ddc
    ddc=DictDBLocal('test_dictdb')

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print('verifying...')
    for a in range(1,250):
        if ddc[a]!=str(a):
            raise RuntimeError

    ddc.standby()

    del ddc

    os.unlink("test_dictdb.dbi")
    print("unlink test")

    ddc=DictDBLocal('test_dictdb')

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print('verifying...')
    for a in range(1,250):
        if ddc[a]!=str(a):
            raise RuntimeError

    ddc.purge()

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print('verifying...')
    for a in range(1,250):
        if ddc[a]!=str(a):
            raise RuntimeError

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    lst.sort()

    for a in lst:
        del ddc[a[1]]

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    del ddc
    ddc=DictDBLocal('test_dictdb')

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    ddc.purge()

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    del ddc
    ddc=DictDBLocal('test_dictdb')

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)

    print(list(ddc.keys()))

    ddc.reset()

    del ddc
    ddc=DictDBLocal('test_dictdb')

    info = ddc.get_info()
    print('used:',info.used,'wasted:',info.wasted,'extant:',info.used-info.wasted)
