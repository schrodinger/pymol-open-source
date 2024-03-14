'''
(c) 2010-2012 Thomas Holder (https://github.com/speleo3/pymol-psico)
(c) 2016 Thomas Holder, Schrodinger, Inc.

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def pdb2pqr(name, selection='all', ff='amber', debump=1, opt=1, assignonly=0,
        ffout='', ph=None, neutraln=0, neutralc=0, state=-1, preserve=0,
        exe='pdb2pqr', quiet=1):
    '''
DESCRIPTION

    Creates a new molecule object from a selection and adds missing atoms,
    assignes charges and radii using PDB2PQR.

    http://www.poissonboltzmann.org/pdb2pqr/

USAGE

    pdb2pqr name [, selection [, ff [, debump [, opt [, assignonly [, ffout [,
        ph [, neutraln [, neutralc [, state [, preserve ]]]]]]]]]]]

ARGUMENTS

    name = string: name of object to create or modify

    selection = string: atoms to include in the new object {default: all}

    ff = string: forcefield {default: amber}
    '''
    debump, opt, assignonly = int(debump), int(opt), int(assignonly)
    neutraln, neutralc = int(neutraln), int(neutralc)
    quiet = int(quiet)

    args = ['--ff=' + ff, '--chain']
    if not debump:
        args.append('--nodebump')
    if not opt:
        args.append('--noopt')
    if assignonly:
        args.append('--assign-only')
    if ffout:
        args.append('--ffout=' + ffout)
    if ph is not None:
        args.append('--with-ph=%f' % float(ph))
    if neutraln:
        args.append('--neutraln')
    if neutralc:
        args.append('--neutralc')
    if not quiet:
        args.append('--verbose')

    r = pdb2pqr_cli(name, selection, args, state, preserve, exe, quiet)

    if not quiet:
        if r:
            print(r)

        print(' pdb2pqr: done')

    return r


def pdb2pqr_cli(name, selection, options, state=-1, preserve=0,
        exe='pdb2pqr', quiet=1, fixrna=0, _proclist=None):
    import os, tempfile, subprocess, shutil

    state, preserve, quiet = int(state), int(preserve), int(quiet)

    if cmd.is_string(options):
        import shlex
        options = shlex.split(options)

    args = [cmd.exp_path(exe)] + list(options)

    tmpdir = tempfile.mkdtemp()
    # Input format should be PDB, but use PQR instead because we can support
    # multi-state assemblies by not writing MODEL records.
    infile = os.path.join(tmpdir, 'in.pqr')
    outfile = os.path.join(tmpdir, 'out.pqr')
    args.extend([infile, outfile])

    # For some reason, catching stdout/stderr with PIPE and communicate()
    # blocks terminate() calls from terminating early. Using a file
    # redirect doesn't show this problem.
    f_stdout = open(os.path.join(tmpdir, 'stdout.txt'), 'w+')

    # RNA resdiue names must be RA, RC, RG, and RU
    tmpmodel = ''
    if int(fixrna) and cmd.count_atoms('(%s) & resn A+C+G+U' % (selection)):
        tmpmodel = cmd.get_unused_name('_tmp')
        cmd.create(tmpmodel, selection, zoom=0)
        cmd.alter(tmpmodel + ' & polymer.nucleic & resn A+C+G+U',
                'resn = "R" + resn')
        selection = tmpmodel

    try:
        cmd.save(infile, selection, state)

        p = subprocess.Popen(args, cwd=tmpdir,
                stdin=subprocess.PIPE,  # Windows pythonw fix
                stdout=f_stdout,
                stderr=f_stdout)
        p.stdin.close()  # Windows pythonw fix

        if _proclist is not None:
            _proclist.append(p)

        p.wait()

        # This allows PyMOL to capture it and display the output in the GUI.
        f_stdout.seek(0)
        print(f_stdout.read().rstrip())

        if p.returncode == -15: # SIGTERM
            raise CmdException('pdb2pqr terminated')

        if p.returncode != 0:
            raise CmdException('%s failed with exit status %d' %
                    (args[0], p.returncode))

        warnings = [L[10:] for L in open(outfile) if L.startswith('REMARK   5')]
        warnings = ''.join(warnings).strip()

        cmd.load(outfile, name)

        return warnings

    except OSError as e:
        print(e)
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        if tmpmodel:
            cmd.delete(tmpmodel)

        f_stdout.close()

        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)


def _is_exe(exe):
    import os, sys
    if os.path.exists(exe):
        return True
    if sys.platform.startswith('win'):
        return os.path.exists(exe + '.exe')
    return False


def prepwizard(name, selection='all', options='', state=-1,
        preserve=0, exe='$SCHRODINGER/utilities/prepwizard', quiet=1,
        _proclist=None):
    '''
DESCRIPTION

    Run the SCHRODINGER Protein Preparation Wizard. Builds missing side
    chains and converts MSE to MET. Other non-default options need to be
    passed with the "options=" argument.

USAGE

    prepwizard name [, selection [, options [, state ]]]

ARGUMENTS

    name = str: name of object to create

    selection = str: atoms to send to the wizard {default: all}

    options = str: additional command line options for prepwizard

    state = int: object state {default: -1 (current)}
    '''
    import os, tempfile, subprocess, shutil, shlex

    state, preserve, quiet = int(state), int(preserve), int(quiet)

    exe = cmd.exp_path(exe)
    if not _is_exe(exe):
        if 'SCHRODINGER' not in os.environ:
            print(' Warning: SCHRODINGER environment variable not set')
        raise CmdException('no such script: ' + exe)

    args = [exe, '-mse', '-fillsidechains', '-WAIT']

    if options:
        if cmd.is_string(options):
            options = shlex.split(options)
        args.extend(options)

    tmpdir = tempfile.mkdtemp()
    infile = 'in.pdb'
    outfile = 'out.mae'
    args.extend([infile, outfile])

    try:
        cmd.save(os.path.join(tmpdir, infile), selection, state)

        env = dict(os.environ)
        env.pop('PYTHONEXECUTABLE', '')  # messes up Python on Mac

        import pymol
        if pymol.IS_WINDOWS:
            # Fix for 2020-4 (PYMOL-3572)
            import ctypes
            ctypes.windll.kernel32.SetDllDirectoryW(None)

        p = subprocess.Popen(args, cwd=tmpdir,
                env=env,
                universal_newlines=True,
                stdin=subprocess.PIPE,  # Windows pythonw fix
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT)

        if _proclist is not None:
            _proclist.append(p)

        print(p.communicate()[0].rstrip())

        if p.wait() == -15: # SIGTERM
            raise CmdException('prepwizard terminated')

        if p.returncode != 0:
            raise CmdException('%s failed with exit status %d' % (args[0], p.returncode))

        cmd.load(os.path.join(tmpdir, outfile), name)
    except OSError as e:
        print(e)
        raise CmdException('Cannot execute "%s"' % (exe))
    finally:
        logfile = os.path.join(tmpdir, 'in.log')
        if os.path.exists(logfile):
            with open(logfile) as handle:
                print(handle.read())
        if not preserve:
            shutil.rmtree(tmpdir)
        elif not quiet:
            print(' Notice: not deleting', tmpdir)

    if not quiet:
        print(' prepwizard: done')


# vi: ts=4:sw=4:smarttab:expandtab
