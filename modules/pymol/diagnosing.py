import sys

def _unicode(s):
    # Return - if possible - unicode, or if the byte string can't be
    # decoded, return an all-ascii representation (which is fine for
    # reporting in diagnostics)
    if isinstance(s, bytes):
        try:
            s = s.decode(sys.getfilesystemencoding() or 'ascii')
        except UnicodeError:
            s = repr(s)
    return s

def diagnostics(filename='', compact=0, quiet=1):
    '''
DESCRIPTION

    Get system level diagnostics data

USAGE

    diagnostics [ filename ]

ARGUMENTS

    filename = str: If given, write output to text file
    '''
    import platform
    import os
    import re
    import glob
    import time
    import textwrap
    from pymol import invocation
    from pymol import cmd, CmdException

    compact, quiet = int(compact), int(quiet)

    TZ = '%+05d' % (time.timezone / 36)

    version = cmd.get_version()
    body = u'PyMOL %s\n' % version[0]

    if not compact:
        if version[3]:
            body += 'build date: %s %s\n' % (time.ctime(version[3]), TZ)
        if version[4]:
            body += 'git sha: %s\n' % version[4]

    try:
        # pymol conda package version
        condameta = glob.glob(os.path.join(sys.prefix, 'conda-meta',
            'pymol-' + version[0] + '*.json'))
        if condameta:
            import json
            d = json.load(open(condameta[0]))
            body += 'conda build: %s %s\n' % (d.get('build'),
                    d.get('schannel') or d.get('channel'))
    except BaseException as e:
        print(e)

    if not compact:
        body += '\nLicense Information:\n'

    body += 'Open-Source Build\n'

    if not compact:
        body += '\nOperating System:\n'
    body += platform.platform() + '\n'
    if platform.system() == 'Linux':
        body += platform.version() + '\n'

    if not compact:
        body += '\nOpenGL Driver:\n'
    renderer = cmd.get_renderer()
    body += str(renderer[0] or '(none)') + '\n'
    body += str(renderer[1] or '(none)') + '\n'
    body += str(renderer[2] or '(none)') + '\n'

    try:
        from pymol.Qt import QtCore
        body += '%s %s (Qt %s)\n' % (
                QtCore.__name__.split('.')[0],
                QtCore.PYQT_VERSION_STR,
                QtCore.QT_VERSION_STR)
    except Exception as e:
        body += '(%s)\n' % (e,)

    if not compact:
        body += '\nPython:\n'
        body += '%s\n' % sys.version
        body += 'prefix=%s\n' % _unicode(sys.prefix)
        body += 'executable=%s\n' % _unicode(sys.executable)
        body += u'filesystemencoding={}\n'.format(sys.getfilesystemencoding())

        body += '\nStartup Scripts:\n'
        pymolrc = invocation.get_user_config()
        if pymolrc:
            pymolrc = map(_unicode, pymolrc)
            body += '\n'.join(pymolrc) + '\n'
        else:
            body += '(no pymolrc file found)\n'

        body += '\nQt, Python and PyMOL Environment Variables:\n'
        for key in sorted(os.environ):
            if re.search(r'^PY|QT|^LANG', key) and not (key == 'PYMOL_SCRIPTS'):
                body += u'%s=%s\n' % (key, _unicode(os.environ[key]))

        body += '\nPATH:\n'
        body += textwrap.fill(_unicode(os.getenv('PATH', '')), 78) + '\n'

        body += '\nDiagnostics collected on %s %s\n' % (time.ctime(), TZ)

    if filename:
        filename = cmd.exp_path(filename)
        if not filename.endswith('.txt'):
            raise CmdException('filename must have .txt extension')

        with open(filename, 'w') as handle:
            handle.write(body)
        print('Diagnostics written to "%s"' % filename)
    elif not quiet:
        print(body.rstrip())

    return body
