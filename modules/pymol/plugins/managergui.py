'''
PyMOL Plugins Engine, Graphical Plugin Manager

(c) 2011-2012 Thomas Holder, PyMOL OS Fellow
License: BSD-2-Clause

'''

import sys
import Pmw
import tkinter as Tkinter
import tkinter.simpledialog as tkSimpleDialog
from .legacysupport import tkMessageBox, tkFileDialog
from . import pref_get

default_pad = {'padx': 5, 'pady': 5}
default_top = {'fill': 'x', 'anchor': 'n', 'side': 'top'}
default_top.update(default_pad)

def confirm_network_access():
    '''
    Popup dialog with network access notification (only once per session)
    '''
    self = confirm_network_access
    if self.ok < 0:
        tkMessageBox.showinfo('Info',
                'Network download has been disabled, sorry!')
        return False
    if self.ok > 0:
        return True
    self.ok = tkMessageBox.askyesno('Confirm', 'PyMOL will now download'
            ' executable code from the internet! Proceed?')
    return self.ok

# valid values: 1=never ask, 0=ask once per session, -1=network access disabled
confirm_network_access.ok = pref_get('network_access_ok', 0)

def manager_dialog():
    '''
    Create the Plugin Manager dialog (with Pmw)
    '''
    from .legacysupport import get_tk_root
    dialog = PluginManager(get_tk_root())

def plugin_info_dialog(parent, info):
    '''
    Shows a popup with plugin info (name, metadata, doc-string, ...)
    '''
    dialog = Pmw.MegaToplevel(parent,
            title = 'Info about plugin ' + info.name)
    grid = dialog.interior()

    highlighted_labels = ['Citation', 'Citation-Required']

    bg_important = '#ff6666'
    bg_notice = '#ffff99'

    def add_line(label, text):
        bg = None
        if label in highlighted_labels:
            bg = bg_important
        row = grid.grid_size()[1]
        Tkinter.Label(grid, text=label + ':', bg=bg).grid(row=row, column=0, sticky='nw', padx=5, pady=2)
        e = Tkinter.Entry(grid)
        e.insert(0, str(text))
        e.config(state='readonly')
        e.grid(row=row, column=1, sticky='nwe', padx=5, pady=2)

    if info.get_citation_required():
        Tkinter.Label(grid, text='This Plugin requires a citation. See below for details.',
                bg=bg_important, padx=10, pady=10).grid(columnspan=2, sticky='nesw')

    add_line('Name', info.name)
    if not info.is_temporary:
        add_line('Python Module Name', info.mod_name)
        add_line('Filename', info.filename)
    grid.columnconfigure(1, weight=1)

    metadata = info.get_metadata()
    for label, value in metadata.items():
        add_line(label, value)

    if not info.is_temporary:
        if not info.loaded:
            Tkinter.Label(grid, text='more information might be available after plugin is loaded',
                    bg=bg_notice, padx=10, pady=10).grid(columnspan=2, sticky='nesw')
        else:
            add_line('commands', ', '.join(info.commands))

    docstring = info.get_docstring()
    if docstring is not None:
        st = Pmw.ScrolledText(grid, text_wrap='none', text_padx = 4, text_pady = 4)
        st.appendtext(docstring.strip())
        st.configure(text_state='disabled')
        st.grid(columnspan=2, sticky='nesw', padx=5, pady=2)
        grid.rowconfigure(grid.grid_size()[1] - 1, weight=1)
    elif info.loaded or info.is_temporary:
        Tkinter.Label(grid, text='no documentation available',
                bg=bg_notice, padx=10, pady=10).grid(columnspan=2, sticky='nesw')

class PluginManager(Pmw.MegaToplevel):
    '''
    Graphical Plugin Manager Dialog
    '''

    def __init__(self, parent):
        self.super = self.__class__.__bases__[0]
        self.super.__init__(self, parent, title='Plugin Manager')

        self.minsize(600, 400)
        master = self.interior()

        # save button (only show if not in "instantsave" mode)
        from . import pref_get, pref_save
        def c_save():
            pref_save()
            b_save.configure(background=b_save_background)
        b_save = Tkinter.Button(master, text='Save Settings', command=c_save)
        b_save_background = b_save.cget('background')
        PluginManager.b_save = b_save
        if not pref_get('instantsave', True):
            b_save.pack(side='bottom', **default_pad)

        # main tabs
        notebook = Pmw.NoteBook(master)
        notebook.pack(fill='both', expand=1, padx=5, pady=5)

        self.page_installed(notebook)
        self.page_install_new(notebook)
        self.page_settings(notebook)
        self.page_about(notebook)

    def page_installed(self, notebook):
        page = notebook.add('Installed Plugins')

        # frames
        f_filter = Tkinter.Frame(page)
        f_installed = InstalledPluginsWidget(page)

        def filter_command():
            f_installed.filter(e_filter.getvalue())

        # filter items
        Tkinter.Checkbutton(f_filter, text='startup', variable=f_installed.v_fstartup, command=filter_command).pack(side='right')
        Tkinter.Checkbutton(f_filter, text='loaded', variable=f_installed.v_floaded, command=filter_command).pack(side='right')
        e_filter = Pmw.EntryField(f_filter, labelpos='w', label_text='Filter:', modifiedcommand=filter_command)
        e_filter.pack(side='left', expand=1, fill='x')

        # enable/disable all items
        f_all = Tkinter.Frame(page)
        Tkinter.Button(f_all, text='startup all', command=f_installed.startup_all).pack(side='left')
        Tkinter.Button(f_all, text='startup none', command=f_installed.startup_none).pack(side='left')

        # pack
        f_filter.pack(**default_top)
        f_all.pack(side='bottom', anchor='w', **default_pad)
        f_installed.pack(side='top', fill='both', expand=1, **default_pad)

        self.f_installed = f_installed

    def page_install_new(self, notebook):
        from pymol import Scratch_Storage, CmdException

        page = notebook.add('Install New Plugin')

        # local file

        def c_local():
            from .legacysupport import installPlugin, get_pmgapp
            installPlugin(get_pmgapp())
            self.f_installed.reload()

        w = Pmw.Group(page, tag_text='Install from local file')
        w.pack(**default_top)

        b_local = Tkinter.Button(w.interior(), text='Choose file...', command=c_local)
        b_local.pack(side='left', **default_pad)

        if confirm_network_access.ok < 0:
            Tkinter.Label(page,
                    text='Network access has been disabled').pack(expand=1)
            return

        # pymolwiki

        def fetchplugin():
            if not confirm_network_access():
                return

            from .installation import installPluginFromFile
            from .repository import fetchscript
            url = e_wiki.get()
            if not len(url):
                return
            import tempfile, shutil
            tmpdir = tempfile.mkdtemp()
            try:
                filename = fetchscript(url, tmpdir, False)
            except CmdException as e:
                tkMessageBox.showerror('Error', 'Fetching Plugin failed.\n' + str(e))
                return

            if filename:
                installPluginFromFile(filename, self.interior())
            shutil.rmtree(tmpdir)
            self.f_installed.reload()

        w = Pmw.Group(page, tag_text='Install from PyMOLWiki or any URL')
        w.pack(**default_top)
        w = w.interior()

        l_repo = Tkinter.Label(w,
                text='Paste a link to a script or plugin, or a PyMOLWiki url which then will be downloaded and '
                'scanned for scripts that extend the PyMOL API', justify='left', anchor='w', wraplength=500)
        l_repo.pack(**default_top)

        e_wiki = Pmw.EntryField(w, labelpos='w', label_text='URL:', command=fetchplugin)
        e_wiki.pack(side='left', fill='x', expand=1, **default_pad)
        b_wiki = Tkinter.Button(w, text='Fetch',
                command=fetchplugin)
        b_wiki.pack(side='right', **default_pad)

        # repository

        w = Pmw.Group(page, tag_text='Install from Repository')
        w.pack(fill='both', **default_pad)

        pw = Pmw.PanedWidget(w.interior(), orient='horizontal')

        pane_left = pw.add('left', min=.2)
        pane_right = pw.add('right', min=.2, max=.5)

        repo_tmp = Scratch_Storage()
        def selecmd_left():
            '''
            Get plugins listing for selected repository.
            '''
            if not confirm_network_access():
                return

            from .repository import guess
            sels = slb_left.getcurselection()
            if len(sels) == 0:
                slb_right.setlist(['- empty -'])
                return
            try:
                url = sels[0]
                repo_tmp.r = guess(url)
                slb_right.setlist(repo_tmp.r.list())
            except:
                slb_right.setlist(['- listing failed -'])

        def infocmd_right():
            '''
            Download file, parse for metadata, show info-popup and delete file
            '''
            from . import PluginInfo
            from .installation import get_name_and_ext, extract_zipfile, zip_extensions
            sels = slb_right.getcurselection()
            if len(sels) == 0:
                return
            import tempfile, shutil, os
            tmpdir = tempfile.mkdtemp()
            tmpdirs = [tmpdir]
            try:
                name = sels[0]
                filename = repo_tmp.r.copy(name, tmpdir)
                name, ext = get_name_and_ext(filename)
                if ext in zip_extensions:
                    tmpdir, dirnames = extract_zipfile(filename, ext)
                    tmpdirs.append(tmpdir)
                    name = dirnames[-1]
                    filename = os.path.join(os.path.join(tmpdir, *dirnames), '__init__.py')
                info = PluginInfo(name, filename)
                plugin_info_dialog(self.interior(), info)
            except:
                tkMessageBox.showinfo('Error', 'Could not get plugin info')
            finally:
                for tmpdir in tmpdirs:
                    shutil.rmtree(tmpdir)

        def selecmd_right():
            '''
            Download plugin from repository and install it.
            '''
            from .installation import installPluginFromFile, get_plugdir
            sels = slb_right.getcurselection()
            if len(sels) == 0:
                return
            plugdir = get_plugdir(self.interior())
            if not plugdir:
                return
            import tempfile, shutil, os
            tmpdir = tempfile.mkdtemp()
            try:
                for name in sels:
                    filename = repo_tmp.r.copy(name, tmpdir)
                    installPluginFromFile(filename, self.interior(), plugdir)
            except:
                err = str(sys.exc_info()[1])
                tkMessageBox.showinfo('Error', 'Could not install plugin ' + name + '\n\n' + err)
            finally:
                shutil.rmtree(tmpdir)
            self.f_installed.reload()

        slb_left = Pmw.ScrolledListBox(pane_left, items=(
                'http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/',
                'https://github.com/Pymol-Scripts/Pymol-script-repo',

                # for testing
                'http://www.thomas-holder.de/projects/pymol/repository/',
                ),
            listbox_height=30, labelpos='nw', label_text='Repositories',
            selectioncommand=selecmd_left)

        slb_right = Pmw.ScrolledListBox(pane_right, items=(),
            listbox_selectmode=Tkinter.EXTENDED,
            listbox_height=30, labelpos='nw', label_text='Items',)

        repo_bb_left = Pmw.ButtonBox(pane_left)
        repo_bb_right = Pmw.ButtonBox(pane_right)

        def dummy_command():
            showinfo = tkMessageBox.showinfo
            showinfo('Dummy', 'Not implemented', parent=self.interior())

        def slb_left_setlist(items):
            slb_left.setlist(items)
            self.b_save.configure(background='red')

        def slb_left_add():
            url = tkSimpleDialog.askstring('Repository URL', 'Please enter Repository URL', parent=self.interior())
            if not url:
                return
            items = list(slb_left.get())
            items.append(url)
            slb_left_setlist(items)
            slb_left.setvalue([url]) # provide list to avoid unicode problem

        def slb_left_remove():
            v = slb_left.getvalue()
            items = list(slb_left.get())
            for item in v:
                items.remove(item)
            slb_left_setlist(items)

        repo_bb_left.add('Add ...', command=slb_left_add)
        repo_bb_left.add('Remove', command=slb_left_remove)
        repo_bb_right.add('Info', command=infocmd_right)
        repo_bb_right.add('Install', command=selecmd_right)

        repo_bb_left.pack(side='bottom', fill='x',)
        repo_bb_right.pack(side='bottom', fill='x',)

        slb_left.pack(fill='both', **default_pad)
        slb_right.pack(fill='both', **default_pad)
        pw.pack(fill='both')

    def page_settings(self, notebook):
        page = notebook.add('Settings')

        from . import get_startup_path, set_startup_path

        # plugin search path
        w = Pmw.Group(page, tag_text='Plugin override search path')
        w.pack(**default_top)

        def slb_path_setlist(items):
            set_startup_path(items)
            slb_path.setlist(items)
            self.b_save.configure(background='red')

        def slb_path_add():
            import os
            from .installation import get_default_user_plugin_path as userpath
            d = tkFileDialog.askdirectory(initialdir=userpath(), parent=self.interior())
            if not len(d):
                return
            if not os.path.exists(d):
                os.makedirs(d)
            items = list(slb_path.get())
            items.append(d)
            slb_path_setlist(items)
            slb_path.setvalue([d]) # provide list to avoid unicode problem

        def slb_path_remove():
            v = slb_path.getvalue()
            items = list(slb_path.get())
            for item in v:
                items.remove(item)
            slb_path_setlist(items)

        def slb_path_move(j=1):
            v = slb_path.getvalue()
            if len(v) == 0:
                return
            items = list(slb_path.get())
            i = items.index(v[0])
            if not (0 <= i+j < len(items)):
                return
            items.pop(i)
            items.insert(i+j, v[0])
            slb_path_setlist(items)
            slb_path.setvalue([v[0]]) # provide list to avoid unicode problem

        slb_path = Pmw.ScrolledListBox(w.interior(),
                items=get_startup_path(True), listbox_height=4)
        if slb_path.get(0):
            slb_path.setvalue([slb_path.get(0)]) # provide list to avoid unicode problem
        slb_path.pack(**default_top)

        bb_path = Pmw.ButtonBox(w.interior())
        bb_path.add('Add new directory...', command=slb_path_add)
        bb_path.add('Remove', command=slb_path_remove)
        bb_path.add('Move up', command=lambda: slb_path_move(-1))
        bb_path.add('Move down', command=slb_path_move)
        bb_path.pack(side='top', fill='x')

        l_path = Tkinter.Label(w.interior(),
                text='PyMOL restart required in order to find plugins in modified plugin search path')
        l_path.pack(**default_top)

        w = Pmw.Group(page, tag_text='Preferences (Read-Only)')
        w.pack(**default_top)
        w = w.interior()

        from . import preferences
        for row, key in enumerate(preferences):
            value = str(preferences[key])
            if len(value) > 100:
                value = value[:97] + '...'
            Tkinter.Label(w, text=key).grid(row=row, column=0, sticky='nw', padx=3, pady=2)
            e = Tkinter.Entry(w)
            e.insert(0, value)
            e.config(state='readonly')
            e.grid(row=row, column=1, sticky='nwe', padx=5, pady=2)
        w.columnconfigure(1, weight=1)

    def page_about(self, notebook):
        page = notebook.add('About')

        st = Pmw.ScrolledText(page, text_wrap='word')
        st.appendtext('Plugins are external modules which extend PyMOL\'s capabilities.\n\n'
                'Plugins can provide new commands and/or add menu items to the "Plugin" menu.\n\n'
                'For technical details, visit\n'
                'http://pymolwiki.org/index.php/PluginArchitecture\n'
                'http://pymolwiki.org/index.php/Script_Tutorial\n')
        st.configure(text_state='disabled')
        st.pack(fill='both', expand=1, **default_pad)

class InstalledPluginsWidget(Pmw.ScrolledFrame):
    '''
    Scrolled widget that shows all installed plugins.
    '''

    def __init__(self, parent):
        self.super = self.__class__.__bases__[0]
        self.super.__init__(self, parent,
                vscrollmode='static',
                horizflex='expand')

        # filter variables
        self.v_floaded  = Tkinter.BooleanVar(parent, False)
        self.v_fstartup = Tkinter.BooleanVar(parent, False)

        self.reload()

    def reload(self):
        from . import plugins

        for child in self.children():
            child.pack_forget()
            child.destroy()

        for info in sorted(plugins.values(), key=lambda i: i.name.lower()):
            PluginWidget(self, info)

        # pack children
        self.filter()

        self.bind_mousewheel()

    def children(self):
        return self.interior().winfo_children()

    def filter(self, text=''):
        '''
        Pack all children that match filter criteria.
        '''
        for child in self.children():
            child.pack_forget()

        for child in self.children():
            if text.lower() not in child.info.name.lower():
                continue
            if self.v_fstartup.get() and not child.info.autoload:
                continue
            if self.v_floaded.get() and not child.info.loaded:
                continue
            child.pack()

    def bind_mousewheel(self):
        '''
        Enable scrolling with mouse wheel.
        '''
        delta = 0.1

        def scrollUp(event=None):
            self.yview('scroll', -delta, 'pages')

        def scrollDown(event=None):
            self.yview('scroll', delta, 'pages')

        def bind_rec(other):
            if hasattr(other, 'interior'):
                other = other.interior()
            other.bind('<Button-4>', scrollUp)
            other.bind('<Button-5>', scrollDown)
            for child in other.winfo_children():
                bind_rec(child)
        bind_rec(self)

    def startup_all(self):
        for child in self.children():
            child.w_startup.select()
            if not child.info.autoload:
                child.c_startup(False)

    def startup_none(self):
        for child in self.children():
            child.w_startup.deselect()
            if child.info.autoload:
                child.c_startup(False)

class PluginWidget(Tkinter.Frame):
    '''
    Item of InstalledPluginsWidget that represents a installed plugin.
    '''

    def __init__(self, parent, info):
        self.info = info

        if hasattr(parent, 'interior'):
            parent = parent.interior()

        self.super = self.__class__.__bases__[0]
        self.super.__init__(self, parent,
                borderwidth=1,
                padx=5, pady=5,
                relief=Tkinter.SUNKEN)

        topmaster = self

        master = Tkinter.Frame(topmaster)
        master.pack(side=Tkinter.TOP, fill='x')

        w_title = Tkinter.Label(master, text=info.name)
        w_title.pack(side=Tkinter.LEFT)

        self.w_version = Tkinter.Label(master, text='', foreground='gray50')
        self.w_version.pack(side=Tkinter.LEFT)

        self.w_loadtime = Tkinter.Label(master, text='', foreground='gray50')
        self.w_loadtime.pack(side=Tkinter.RIGHT)

        master = Tkinter.Frame(topmaster)
        master.pack(side=Tkinter.TOP, fill='x')

        self.w_info = Tkinter.Button(master, text='Info', command=self.plugin_info)
        self.w_info.pack(side=Tkinter.LEFT)

        self.w_enable = Tkinter.Button(master, text='Load', command=self.plugin_load)
        self.w_enable.pack(side=Tkinter.LEFT)

        self.v_startup = Tkinter.BooleanVar(master, self.info.autoload)
        self.w_startup = Tkinter.Checkbutton(master, text='Load on startup',
                variable=self.v_startup, command=self.c_startup)
        self.w_startup.pack(side=Tkinter.LEFT)

        self.f_settings = Tkinter.Frame(master)
        self.f_settings.pack(side=Tkinter.LEFT)

        w_uninstall = Tkinter.Button(master, text='Uninstall', command=self.plugin_remove)
        w_uninstall.pack(side=Tkinter.RIGHT)

        self.status_update()

    def c_startup(self, askload=True):
        '''
        Callback for "autoload" checkbox.
        '''
        self.info.autoload = self.v_startup.get()
        if askload and self.info.autoload and not self.info.loaded:
            if tkMessageBox.askyesno('Confirm', 'Load plugin now?', parent=self):
                self.plugin_load()
        PluginManager.b_save.configure(background='red')

    def pack(self):
        self.super.pack(self, fill='x', anchor=Tkinter.N)

    def plugin_info(self):
        plugin_info_dialog(self, self.info)

    def plugin_load(self):
        self.info.load()
        self.status_update()

    def plugin_remove(self):
        if self.info.uninstall(self):
            self.pack_forget()
            self.destroy()

    def status_update(self):
        '''
        Update contents of this item that may have changed after loading plugin.
        '''
        version = self.info.get_version()
        if version:
            text = 'version ' + self.info.get_version()
            self.w_version.config(text=text)

        if self.info.loaded:
            text = 'Took %.3f seconds to load' % (self.info.loadtime)
            self.w_loadtime.config(text=text)
            self.w_enable.config(state=Tkinter.DISABLED)

            # add 'Settings' button
            if hasattr(self.info.module, 'settings_dialog') \
                    and not len(self.f_settings.winfo_children()):
                Tkinter.Button(self.f_settings, text='Settings',
                        command=self.info.module.settings_dialog).pack()

# vi:expandtab:smarttab:sw=4:nowrap
