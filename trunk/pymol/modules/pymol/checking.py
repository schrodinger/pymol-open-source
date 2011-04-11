
import types
import cmd

def _raising(code=-1,_self=cmd):
    # WARNING: internal routine, subject to change
    if isinstance(code, types.IntType):
        if code<0:
            return _self.get_setting_legacy("raise_exceptions")
    return 0

def is_string(obj):
    return (isinstance(obj,types.StringType) or isinstance(obj,types.UnicodeType))

def is_list(obj):
    return isinstance(obj,types.ListType)

def is_dict(obj):
    return isinstance(obj,types.DictType)

def is_tuple(obj):
    return isinstance(obj,types.TupleType)

def is_sequence(obj):
    return isinstance(obj,types.ListType) or isinstance(obj,types.TupleType)

def is_error(result): # errors are always negative numbers
    if isinstance(result,types.IntType):
        return (result<0)
    return 0

def is_ok(result): # something other than a negative number
    if isinstance(result,types.IntType):
        return (result>=0)
    return 1
