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

from chempy import Storage

if True:
    import pickle

    # Python 3: Unpickle printable ASCII strings in [TAB, DEL) to unicode,
    # and everything else to bytes.
    # This uses the Python implementation of pickle, not the fast cpython
    # version, which unfortunately is much slower.

    def _decode_string(self, value):
        if all(8 < b < 127 for b in value):
            return value.decode('ascii')
        return value

    pickle._Unpickler._decode_string = _decode_string

    # string pickling backported from Python 2

    def save_str(self, obj):
        from struct import pack

        if self.bin:
            if not isinstance(obj, bytes):
                obj = obj.encode('utf-8', 'ignore')
            n = len(obj)
            if n < 256:
                self.write(pickle.SHORT_BINSTRING + pack("<B", n) + obj)
            else:
                self.write(pickle.BINSTRING + pack("<i", n) + obj)
        else:
            obj_repr = repr(obj).encode('utf-8', 'ignore').lstrip(b'b')
            self.write(pickle.STRING + obj_repr + b'\n')

    pickle._Pickler.dispatch[str] = save_str
    pickle._Pickler.dispatch[bytes] = save_str

    class cPickle:
        dumps = pickle.dumps
        dump = pickle.dump
        load = pickle._load

        def loads(s):
            if not isinstance(s, bytes):
                s = s.encode(errors='ignore')
            return pickle._loads(s)

        @classmethod
        def configure_legacy_dump(cls, py2=False):
            if py2:
                cls.dumps = pickle._dumps
                cls.dump = pickle._dump
            else:
                cls.dumps = pickle.dumps
                cls.dump = pickle.dump


class PKL(Storage):

    def fromFile(self,fname,**params):
        fp = self.my_open(fname,'rb')
        result = cPickle.load(fp)
        fp.close()
        return result

#---------------------------------------------------------------------------
    def toFile(self,indexed,fname,**params):
        fp = open(fname,'wb')
        if 'bin' not in params:
            result = cPickle.dump(indexed,fp,1)
        else:
            result = cPickle.dump(indexed,fp,params['bin'])
        fp.close()

#---------------------------------------------------------------------------
    def fromStream(self,fp,**params):
        try:
            return cPickle.load(fp)
        except EOFError:
            return None

#---------------------------------------------------------------------------
    def toStream(self,indexed,fp,**params):
        if 'bin' not in params:
            result = cPickle.dump(indexed,fp,1)
        else:
            result = cPickle.dump(indexed,fp,params['bin'])

#---------------------------------------------------------------------------
    def fromString(self,st):
        return cPickle.loads(st)

#---------------------------------------------------------------------------
    def toString(self,model):
        return cPickle.dumps(model)
