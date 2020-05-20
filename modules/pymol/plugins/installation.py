'''
PyMOL Plugins Engine, Installation Routines

(c) 2011-2012 Thomas Holder, PyMOL OS Fellow
License: BSD-2-Clause

'''

import os

# supported file types for installation. Do not support pyc and pyo binaries,
# we want text files that can be parsed for metadata.
zip_extensions = ['zip', 'tar.gz']
supported_extensions = ['py'] + zip_extensions

class InstallationCancelled(Exception):
    pass

class BadInstallationFile(Exception):
    pass

def get_default_user_plugin_path():
    '''
    User plugin directory defaults to ~/.pymol/startup on Linux and to
    %APPDATA%\pymol\startup on windows.
    '''
    if 'APPDATA' in os.environ:
        return os.path.join(os.environ['APPDATA'], 'pymol', 'startup')
    return os.path.expanduser('~/.pymol/startup')

def is_writable(dirname):
    '''
    Return True if directory is writable.
    '''
    path = os.path.join(dirname, '__check_writable')
    try:
        f = open(path, 'wb')
        f.close()
        os.remove(path)
        return True
    except (IOError, OSError):
        return False

def cmp_version(v1, v2):
    '''
    Compares two version strings. An empty version string is always considered
    smaller than a non-empty version string.

    Uses distutils.version.StrictVersion to evaluate non-empty version strings.
    '''
    if v1 == v2:
        return 0
    if v1 == '':
        return -1
    if v2 == '':
        return 1
    try:
        from distutils.version import StrictVersion as Version
        return cmp(Version(v1), Version(v2))
    except:
        print(' Warning: Version parsing failed for', v1, 'and/or', v2)
        return 0

def get_name_and_ext(ofile):
    '''
    Given a filename, return module name and file extension.

    Examples:
    foo-1.0.py -> ('foo', 'py')
    /foo/bar.tar.gz -> ('bar', 'tar.gz')
    '''
    import re

    basename = os.path.basename(ofile)
    pattern = r'(\w+).*\.(%s)$' % '|'.join(supported_extensions)
    m = re.match(pattern, basename, re.IGNORECASE)

    if m is None:
        raise BadInstallationFile('Not a valid plugin filename (%s).' % (basename))

    return m.group(1), m.group(2).lower()

def check_valid_name(name):
    '''
    Check if "name" is a valid python module name.
    '''
    if '.' in name:
        raise BadInstallationFile('name must not contain dots (%s).' % repr(name))

def extract_zipfile(ofile, ext):
    '''
    Extract zip file to temporary directory
    '''
    if ext == 'zip':
        import zipfile
        zf = zipfile.ZipFile(ofile)
    else:
        import tarfile
        zf = tarfile.open(ofile)
        zf.namelist = zf.getnames
    # make sure pathnames are not absolute
    cwd = os.getcwd()
    namelist = zf.namelist()
    for f in namelist:
        f = os.path.normpath(f)
        if not os.path.abspath(f).startswith(cwd):
            raise BadInstallationFile('ZIP file contains absolute path names')
    # analyse structure
    namedict = dict()
    for f in namelist:
        x = namedict
        for part in f.split('/'): # even on windows this is a forward slash (not os.sep)
            if part != '':
                x = x.setdefault(part, {})
    if len(namedict) == 0:
        raise BadInstallationFile('Archive empty.')

    # case 1: zip/<name>/__init__.py
    names = [(name,)
            for name in namedict
            if '__init__.py' in namedict[name]]
    if len(names) == 0:
        # case 2: zip/<name>-<version>/<name>/__init__.py
        names = [(pname, name)
                for (pname, pdict) in namedict.items()
                for name in pdict
                if '__init__.py' in pdict[name]]

    if len(names) == 0:
        raise BadInstallationFile('Missing __init__.py')
    if len(names) > 1:
        # filter out "tests" directory
        names = [n for n in names if n[-1] != 'tests']
    if len(names) > 1:
        raise BadInstallationFile('Archive must contain a single package.')
    check_valid_name(names[0][-1])

    # extract
    import tempfile
    tempdir = tempfile.mkdtemp()
    zf.extractall(tempdir)

    return tempdir, names[0]

def get_plugdir(parent=None):
    '''
    Get plugin directory, ask user if startup path has more than one entry
    '''
    from . import get_startup_path
    plugdirs = get_startup_path()

    if len(plugdirs) == 1:
        return plugdirs[0]

    import sys
    if 'pmg_qt.mimic_tk' in sys.modules:
        from pymol.Qt import QtWidgets
        value, result = QtWidgets.QInputDialog.getItem(None,
            'Select plugin directory',
            'In which directory should the plugin be installed?', plugdirs)
        return value if result else ''

    dialog_selection = []
    def plugdir_callback(result):
        if result == 'OK':
            dialog_selection[:] = dialog.getcurselection()
        dialog.destroy()

    import Pmw
    dialog = Pmw.SelectionDialog(parent, title='Select plugin directory',
            buttons = ('OK', 'Cancel'), defaultbutton='OK',
            scrolledlist_labelpos='n',
            label_text='In which directory should the plugin be installed?',
            scrolledlist_items=plugdirs,
            command=plugdir_callback)
    dialog.component('scrolledlist').selection_set(0)

    # wait for dialog to be closed
    dialog.wait_window()

    if not dialog_selection:
        return ''

    return dialog_selection[0]

def installPluginFromFile(ofile, parent=None, plugdir=None):
    '''
    Install plugin from file.

    Takes python (.py) files and archives which contain a python module.
    '''
    import shutil
    from . import startup, PluginInfo
    from . import get_startup_path, set_startup_path, pref_get
    from .legacysupport import tkMessageBox, get_tk_focused

    if parent is None:
        parent = get_tk_focused()

    showinfo = tkMessageBox.showinfo
    askyesno = tkMessageBox.askyesno

    plugdirs = get_startup_path()

    if not plugdir:
        plugdir = get_plugdir()
        if not plugdir:
            return

    if not is_writable(plugdir):
        user_plugdir = get_default_user_plugin_path()
        if not askyesno('Warning',
                'Unable to write to the plugin directory.\n'
                'Should a user plugin directory be created at\n' + user_plugdir + '?',
                parent=parent):
            showinfo('Error', 'Installation aborted', parent=parent)
            return

        if not os.path.exists(user_plugdir):
            try:
                os.makedirs(user_plugdir)
            except OSError:
                showinfo('Error', 'Could not create user plugin directory', parent=parent)
                return

        plugdir = user_plugdir

    if plugdir not in plugdirs:
        set_startup_path([plugdir] + get_startup_path(True))

    def remove_if_exists(pathname, ask):
        '''
        Remove existing plugin files before reinstallation. Will not remove
        files if installing into different startup directory.
        '''
        if not os.path.exists(pathname):
            return

        is_dir = os.path.isdir(pathname)

        if ask:
            if is_dir:
                msg = 'Directory "%s" already exists, overwrite?' % pathname
            else:
                msg = 'File "%s" already exists, overwrite?' % pathname
            if not tkMessageBox.askyesno('Confirm', msg, parent=parent):
                raise InstallationCancelled('will not overwrite "%s"' % pathname)

        if is_dir:
            shutil.rmtree(pathname)
        else:
            os.remove(pathname)

    def check_reinstall(name, pathname):
        from . import plugins

        if name not in plugins:
            remove_if_exists(pathname, True)
            return

        v_installed = plugins[name].get_version()
        v_new = PluginInfo(name, ofile).get_version()
        c = cmp_version(v_new, v_installed)
        if c > 0:
            msg = 'An older version (%s) of this plugin is already installed. Install version %s now?' % (v_installed, v_new)
        elif c == 0:
            msg = 'Plugin already installed. Reinstall?'
        else:
            msg = 'A newer version (%s) of this plugin is already installed. Install anyway?' % (v_installed)

        if not tkMessageBox.askokcancel('Confirm', msg, parent=parent):
            raise InstallationCancelled

        remove_if_exists(pathname, False)

    name = "unknown" # fallback for error message

    temppathnames = []
    try:
        name, ext = get_name_and_ext(ofile)

        if ext in zip_extensions:
            # import archive

            tempdir, dirnames = extract_zipfile(ofile, ext)
            temppathnames.append((tempdir, 1))

            # install
            name = dirnames[-1]
            odir = os.path.join(tempdir, *dirnames)
            ofile = os.path.join(odir, '__init__.py')
            mod_dir = os.path.join(plugdir, name)
            check_reinstall(name, mod_dir)
            check_valid_name(name)
            shutil.copytree(odir, mod_dir)

            mod_file = os.path.join(mod_dir, '__init__.py')

        elif name == '__init__':
            # import directory
            odir = os.path.dirname(ofile)
            name = os.path.basename(odir)
            mod_dir = os.path.join(plugdir, name)
            check_reinstall(name, mod_dir)
            check_valid_name(name)
            shutil.copytree(odir, mod_dir)

            mod_file = os.path.join(mod_dir, '__init__.py')

        elif ext == 'py':
            # import python file
            mod_file = os.path.join(plugdir, name + '.py')
            check_reinstall(name, mod_file)
            check_valid_name(name)
            shutil.copy(ofile, mod_file)

        else:
            raise UserWarning('this should never happen')

    except InstallationCancelled:
        showinfo('Info', 'Installation cancelled', parent=parent)
        return

    except Exception as e:
        if pref_get('verbose', False):
            import traceback
            traceback.print_exc()
        msg = 'Unable to install plugin "{}".\n{}'.format(name, e)
        showinfo('Error', msg, parent=parent)
        return

    finally:
        for (pathname, is_dir) in temppathnames:
            if is_dir:
                shutil.rmtree(pathname)
            else:
                os.remove(pathname)

    prefix = startup.__name__
    info = PluginInfo(name, mod_file, prefix + '.' + name)

    if info.load(force=1):
        showinfo('Success', 'Plugin "%s" has been installed.' % name, parent=parent)
    else:
        showinfo('Error', 'Plugin "%s" has been installed but initialization failed.' % name, parent=parent)

    if info.get_citation_required():
        if askyesno('Citation Required', 'This plugin requires citation. Show information now?'
                '\n\n(You can always get this information from the Plugin Manager, click the "Info" button there)',
                parent=parent):
            from .managergui import plugin_info_dialog
            plugin_info_dialog(parent, info)

# vi:expandtab:smarttab:sw=4
