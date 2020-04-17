#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC.
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

if True:

    from . import selector
    from .cmd import _cmd,lock,unlock,Shortcut,QuietException, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    cmd = __import__("sys").modules["pymol.cmd"]
    import threading
    import pymol
    import string

    def get_bond_print(obj,max_bond,max_type,_self=cmd):
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.get_bond_print(_self._COb,str(obj),int(max_bond),int(max_type))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def spheroid(object="",average=0,_self=cmd):  # EXPERIMENTAL
        '''
DESCRIPTION

    "spheroid" averages trajectory frames together to create
    an ellipsoid-like approximation of the actual anisotropic
    motion exhibited by the atom over a series of trajectory frames.

USAGE

    spheroid object,average

    average = number of states to average for each resulting spheroid state

    '''
        print("Warning: 'spheroid' is experimental, incomplete, and unstable.")
        with _self.lockcm:
            r = _cmd.spheroid(_self._COb,str(object),int(average))
        return r

    def mem(_self=cmd):
        '''
DESCRIPTION

    "mem" Dumps current memory state to standard output. This is a
    debugging feature, not an official part of the API.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.mem(_self._COb)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def check(selection=None, preserve=0):
        '''
DESCRIPTION

    "check" is unsupported command that may eventually have something
    to do with assigning forcefield parameters to a selection of
    atoms.
    
'''
        # This function relies on code that is not currently part of PyMOL/ChemPy
        # NOTE: the realtime module relies on code that is not yet part of PyMOL/ChemPy
        from chempy.tinker import realtime
        if selection is None:
            arg = cmd.get_names("objects")
            arg = arg[0:1]
            if arg:
                if len(arg):
                    selection = arg
        if selection is not None:
            selection = selector.process(selection)
            realtime.assign("("+selection+")",int(preserve))
            realtime.setup("("+selection+")")

    def fast_minimize(*args, **kwargs):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with doing a quick clean up of the
    molecular structure.
    
'''
        kwargs['_setup'] = 0
        return minimize(*args, **kwargs)

    def minimize(sele='', iter=500, grad=0.01, interval=50, _setup=1, _self=cmd):
        '''
DESCRIPTION

    "fast_minimize" is an unsupported nonfunctional command that may
    eventually have something to do with minimization.
    
'''
        from chempy.tinker import realtime

        if not sele:
            names = _self.get_names("objects")
            if not names:
                return
            sele = names[0]
        sele = '(' + sele + ')'

        if not int(_setup) or realtime.setup(sele):
            _self.async_(realtime.mini, int(iter), float(grad), int(interval), sele)
        else:
            print(" minimize: missing parameters, can't continue")


    def dump(fnam, obj, state=1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    The dump command writes the geometry of an isosurface, isomesh,
    isodot, or map object to a simple text file. Each line contains one
    vertex in case of representations, or one grid point in case of a map.

    For surface objects, XYZ coordinates and the normal are exported.
    Three lines make one triangle (like GL_TRIANGLES).

    For mesh objects, XYZ coordinates are exported (no normals).
    The vertices form line strips (like GL_LINE_STRIP), a blank
    line starts a new strip.

    For dot objects, XYZ coordinates are exported.

    For map objects, XYZ coordinates and the value at the point are
    exported. This forms a grid map.

USAGE

    dump filename, object, state=1, quiet=1

ARGUMENTS

    filename = str: file that will be written
    object = str: object name

EXAMPLE

    fetch 1ubq, mymap, type=2fofc, async=0

    dump gridmap.txt, mymap

    isosurface mysurface, mymap
    dump surfacegeometry.txt, mysurface

    isomesh mymesh, mymap
    dump meshgeometry.txt, mymesh

    isodot mydot, mymap, quiet=1
    dump dotgeometry.txt, mydot

SEE ALSO

    COLLADA export

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.dump(_self._COb, str(fnam), obj, int(state) - 1, int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r


    def dummy(*arg):
        return None

    def test(group=0,index=0,_self=cmd): # generic test routine for development
        '''
DESCRIPTION

    "dump" is an unsupported internal command.

    '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r=_cmd.test(_self._COb,int(group),int(index))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def load_coords(model, oname, state=1): # UNSUPPORTED
        '''
        WARNING: buggy argument list, state get's decremented twice!
        '''
        return pymol.importing.load_coordset(model, oname, int(state)-1)

    def focal_blur(aperture=2.0, samples=10, ray=0, filename='', quiet=1, _self=cmd):
        '''
DESCRIPTION

    Creates fancy figures by introducing a focal blur to the image.
    The object at the origin will be in focus.

USAGE

    focal_blur [ aperture [, samples [, ray [, filename ]]]]

ARGUMENTS

    aperture = float: aperture angle in degrees {default: 2.0}

    samples = int: number of images for averaging {default: 10}

    ray = 0/1: {default: 0}

    filename = str: write image to file {default: temporary}

AUTHORS

    Jarl Underhaug, Jason Vertrees and Thomas Holder

EXAMPLES

    focal_blur 3.0, 50
        '''
        raise pymol.IncentiveOnlyException()

    def callout(name, label, pos='', screen='auto', state=-1, color='front',
            quiet=1, _self=cmd):
        '''
DESCRIPTION

    Create a new screen-stabilized callout object.

ARGUMENTS

    name = str: object name

    label = str: label text

    pos = str or list: anchor in model space as 3-float coord list or atom
    selection. If empty, don't draw an arrow. {default: }

    screen = str or list: position on screen as 2-float list between [-1,-1]
    (lower left) and [1,1] (upper right) or "auto" for smart placement.
    {default: auto}
        '''
        raise pymol.IncentiveOnlyException()

    def desaturate(selection="all", a=0.5, quiet=1, _self=cmd):
        '''
DESCRIPTION

    Desaturate the colors in the given selection.

ARGUMENTS

    selection = str: atom selection {default: all}

    a = float [0..1]: desaturation factor {default: 0.5}
        '''
        raise pymol.IncentiveOnlyException()
