'''
APBS wrapper

(c) 2012 Thomas Holder (https://github.com/speleo3/pymol-psico)
(c) 2016 Thomas Holder, Schrodinger, Inc.

License: BSD-2-Clause
'''

import os

from pymol import cmd, CmdException

template_apbs_in = '''
read
    mol pqr "{pqrfile}"
end
elec
    mg-auto
    mol 1

    fgcent {fgcent} # fine grid center
    cgcent mol 1    # coarse grid center
    fglen {fglen}
    cglen {cglen}
    dime {dime}
    lpbe          # l=linear, n=non-linear Poisson-Boltzmann equation
    bcfl sdh      # "Single Debye-Hueckel" boundary condition
    pdie 2.0      # protein dielectric
    sdie 78.0     # solvent dielectric
    chgm spl2     # Cubic B-spline discretization of point charges on grid
    srfm smol     # smoothed surface for dielectric and ion-accessibility coefficients
    swin 0.3
    temp 310.0    # temperature
    sdens 10.0
    calcenergy no
    calcforce no
    srad {srad}   # solvent radius

    ion charge +1 conc 0.15 radius 2.0
    ion charge -1 conc 0.15 radius 1.8

    write pot dx "{mapfile}"
end
quit
'''


def find_apbs_exe():
    import shutil
    exe = shutil.which('apbs')
    if not exe:
        exe = cmd.exp_path('$SCHRODINGER/utilities/apbs')
        if not os.path.exists(exe):
            return None
    return exe


def validate_apbs_exe(exe):
    '''Get and validate apbs executable.
    Raise CmdException if not found or broken.'''
    import subprocess

    if exe:
        exe = cmd.exp_path(exe)
    else:
        exe = find_apbs_exe() or 'apbs'

    try:
        r = subprocess.call([exe, "--version"],
                stdin=subprocess.PIPE,  # Windows pythonw fix
                stdout=open(os.devnull, "w"), stderr=subprocess.STDOUT)
        if r < 0:
            raise CmdException("Broken executable: " + exe)
    except OSError as e:
        print(e)
        raise CmdException("Cannot execute: " + exe)

    return exe

def map_new_apbs(name, selection='all', grid=0.5, buffer=10.0, state=1,
        preserve=0, exe='', assign=0, focus='', quiet=1, _template='',
        _proclist=None):
    '''
DESCRIPTION

    Create electrostatic potential map with APBS.

    "selection" needs partial charges (partial_charge) and radii (elec_radius)
    assigned. These can be loaded for example from a PQR file.

SEE ALSO

    map_new (coulomb), psico.electrostatics
    '''
    import tempfile, shutil, glob, subprocess

    selection = '(%s) and not solvent' % (selection)
    grid, buffer, state = float(grid), float(buffer), int(state)
    preserve, assign, quiet = int(preserve), int(assign), int(quiet)

    # path to apbs executable
    exe = validate_apbs_exe(exe)

    # temporary directory
    tempdir = tempfile.mkdtemp()
    if not quiet:
        print(' Tempdir:', tempdir)

    # filenames
    pqrfile = os.path.join(tempdir, 'mol.pqr')
    infile = os.path.join(tempdir, 'apbs.in')
    stem = os.path.join(tempdir, 'map')

    # temporary selection
    tmpname = cmd.get_unused_name('_sele')
    cmd.select(tmpname, selection, 0)

    cmd.save(pqrfile, tmpname, state, format='pqr', quiet=quiet)

    # grid dimensions
    extent = cmd.get_extent(tmpname, state)
    extentfocus = cmd.get_extent(focus) if focus else extent
    fglen = [(emax-emin + 2*buffer) for (emin, emax) in zip(*extentfocus)]
    cglen = [(emax-emin + 4*buffer) for (emin, emax) in zip(*extent)]
    # make coarse grid a cube (better for non-globular shapes)
    cglen = [max(cglen)] * 3

    cmd.delete(tmpname)

    apbs_in = {
        'pqrfile': pqrfile,
        'fgcent': 'mol 1',
        'fglen': '%f %f %f' % tuple(fglen),
        'cglen': '%f %f %f' % tuple(cglen),
        'srad': cmd.get('solvent_radius'),
        'mapfile': stem,
    }

    if focus:
        apbs_in['fgcent'] = '%f %f %f' % tuple((emax + emin) / 2.0
                for (emin, emax) in zip(*extentfocus))

    try:
        # apbs will fail if grid does not fit into memory
        # -> on failure repeat with larger grid spacing
        for _ in range(3):
            dime = [1 + max(64, n / grid) for n in fglen]
            apbs_in['dime'] = '%d %d %d' % tuple(dime)

            # apbs input file
            with open(infile, 'w') as f:
                f.write((_template or template_apbs_in).format(**apbs_in))

            # run apbs
            p = subprocess.Popen([exe, infile], cwd=tempdir)

            if _proclist is not None:
                _proclist.append(p)

            r = p.wait()

            if r == -15: # SIGTERM
                raise CmdException('apbs terminated')

            if r == 0:
                break

            if r in (-6, -9):
                grid *= 2.0
                if not quiet:
                    print(' Warning: retry with grid =', grid)
                continue

            raise CmdException('apbs failed with code ' + str(r))

        dx_list = glob.glob(stem + '*.dx')
        if not dx_list:
            dx_list = glob.glob(stem + '*.dxbin')
        if len(dx_list) != 1:
            raise CmdException('dx file missing')

        # load map
        cmd.load(dx_list[0], name, quiet=quiet)
    except OSError as e:
        print(e)
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if not preserve:
            shutil.rmtree(tempdir)
        elif not quiet:
            print(' Notice: not deleting %s' % (tempdir))


# vi: ts=4:sw=4:smarttab:expandtab
