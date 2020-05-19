'''
Experimental MMTF (Macromolecular Transmission Format) I/O library
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sys
import itertools
import numpy
import struct

try:
    import msgpack
except ImportError:
    import umsgpack as msgpack

if True:
    from urllib.request import urlopen
    izip = zip
    izip_longest = itertools.zip_longest
    as_msgpack_key = lambda k: k if isinstance(k, bytes) else k.encode()
    buffer = lambda s, i=0: memoryview(s)[i:]

# should be replaced with a more efficient numpy-array aware iterator
def simpleiter(iterable):
    if isinstance(iterable, numpy.ndarray):
        # Iterating over a numpy array is unreasonably slow. Iterating
        # over it's list copy is much faster!
        return iter(iterable.tolist())
    return iter(iterable)

def asarray(arr, dtype='i'):
    if hasattr(arr, '__len__'):
        return numpy.asarray(arr, dtype)
    return numpy.fromiter(arr, dtype)

MMTF_ENDIAN = '>' # big-endian

########### ENCODINGS ################

class RunLength:
    @staticmethod
    def encode(iterable):
        in_iter = simpleiter(iterable)
        curr = next(in_iter)
        counter = 1
        for item in in_iter:
            if item == curr:
                counter += 1
            else:
                yield curr
                yield counter
                curr = item
                counter = 1
        yield curr
        yield counter

    @staticmethod
    def decode(iterable):
        in_iter = simpleiter(iterable)
        out = []
        extend = out.extend
        for item in in_iter:
            extend([item] * next(in_iter))
        return out

class Delta:
    @staticmethod
    def encode(iterable):
        return numpy.diff(asarray(iterable, 'i4'))

    @staticmethod
    def decode(iterable):
        return asarray(iterable).cumsum(dtype='i4')

class RecursiveIndex:
    def __init__(self, min, max):
        self.limits = (min, max)

    def encode(self, iterable):
        min, max = self.limits
        for curr in simpleiter(iterable):
            while curr >= max:
                yield max
                curr -=  max
            while curr <= min:
                yield min
                curr -= min
            yield curr

    def decode(self, iterable):
        min, max = self.limits
        decoded_val = 0
        for item in simpleiter(iterable):
            decoded_val += item
            if item != max and item != min:
                yield decoded_val
                decoded_val = 0

class IntegerFloats:
    def __init__(self, factor):
        self.factor = factor

    def encode(self, in_floats):
        return (asarray(in_floats, 'f4') * self.factor).astype('i4')

    def decode(self, in_ints):
        return asarray(in_ints, 'f4') / self.factor

class IntegerChars:
    @staticmethod
    def encode(in_chars):
        return [ord(x) for x in in_chars]

    @staticmethod
    def decode(in_ints):
        return [(chr(x) if x else '') for x in simpleiter(in_ints)]

######## BUFFERS ###########

class NumbersBuffer:
    def __init__(self, basetype='i', dectype=''):
        self.enctype = numpy.dtype(MMTF_ENDIAN + basetype)
        self.dectype = numpy.dtype(dectype or basetype)

    def decode(self, in_bytes):
        return numpy.frombuffer(in_bytes, self.enctype).astype(self.dectype)

    def encode(self, in_ints):
        return asarray(in_ints, self.enctype).tostring()

class StringsBuffer:
    def __init__(self, nbytes, encoding='ascii'):
        self.enctype = numpy.dtype('S' + str(nbytes))
        self.encoding = encoding

    def decode(self, in_bytes):
        bstrings = numpy.frombuffer(in_bytes, self.enctype)
        return [b.decode(self.encoding) for b in bstrings]

    def encode(self, strings):
        bstrings = numpy.fromiter((s.encode(self.encoding) for s in strings),
                self.enctype, len(strings))
        return bstrings.tostring()

########## STRATEGIES #############

def _PackedIntBufStrategy(nbytes=1, dectype='i4'):
    m = 1 << (nbytes * 8 - 1)
    return [
        NumbersBuffer('i' + str(nbytes), dectype),
        RecursiveIndex(-m, m - 1),
    ]

strategies = {
    1: [NumbersBuffer('f4')],
    2: [NumbersBuffer('i1')],
    3: [NumbersBuffer('i2')],
    4: [NumbersBuffer('i4')],
    5: lambda length: [StringsBuffer(length)],
    6: [NumbersBuffer('i4'), RunLength, IntegerChars],
    7: [NumbersBuffer('i4'), RunLength],
    8: [NumbersBuffer('i4'), RunLength, Delta],
    9: lambda factor: [NumbersBuffer('i4'), RunLength, IntegerFloats(factor)],
   10: lambda factor: _PackedIntBufStrategy(2) + [Delta, IntegerFloats(factor)],
   11: lambda factor: [NumbersBuffer('i2'), IntegerFloats(factor)],
   12: lambda factor: _PackedIntBufStrategy(2) + [IntegerFloats(factor)],
   13: lambda factor: _PackedIntBufStrategy(1) + [IntegerFloats(factor)],
   14: _PackedIntBufStrategy(2),
   15: _PackedIntBufStrategy(1),
}

# optional parameters format (defaults to 'i' -> one int32 argument)
strategyparamsfmt = {
}

########## MEDIUM LEVEL ARRAY ENCODE/DECODE API ##############

def encode(arr, codec, param=0):
    strategy = strategies[codec]

    if not isinstance(strategy, list):
        strategy = strategy(param)

    buf = struct.pack(MMTF_ENDIAN + 'iii', codec, len(arr), param)

    for handler in reversed(strategy):
        arr = handler.encode(arr)

    buf += arr
    return buf

def encode_int(arr):
    '''Find the best compression for a 32bit int array'''
    return min((encode(arr, codec) for codec in (4, 7, 8)), key=len)

def decode(value):
    codec, length = struct.unpack(MMTF_ENDIAN + 'ii', value[:8])

    strategy = strategies[codec]
    if not isinstance(strategy, list):
        fmt = strategyparamsfmt.get(codec, 'i')
        params = struct.unpack(MMTF_ENDIAN + fmt, value[8:12])
        strategy = strategy(*params)

    value = buffer(value, 12)
    for handler in strategy:
        value = handler.decode(value)

    return value

############### HIGH LEVEL READER API ################

class MmtfReader:
    def __init__(self, data):
        if isinstance(data, bytes):
            if data[:2] != b'\x1f\x8b': # gzip magic number
                self._data = msgpack.unpackb(data)
                return

            import io, gzip
            data = gzip.GzipFile(fileobj=io.BytesIO(data))

        self._data = msgpack.unpack(data)

    @classmethod
    def from_url(cls, url):
        handle = open(url, 'rb') if os.path.isfile(url) else urlopen(url)
        return cls(handle.read())

    def get(self, key, default=None):
        key = as_msgpack_key(key)
        try:
            value = self._data[key]
        except KeyError:
            return default

        if not (key.endswith(b'List') and isinstance(value, bytes)):
            return value

        return decode(value)

    def get_iter(self, key, default=()):
        return simpleiter(self.get(key, default))

    def get_table_iter(self, keys, defaults=None):
        if defaults is None:
            return izip_longest(*[self.get_iter(k) for k in keys])
        return izip(*[self.get_iter(k, itertools.repeat(d))
            for (k, d) in zip(keys, defaults)])
