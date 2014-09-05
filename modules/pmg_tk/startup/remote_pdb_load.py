'''
PDB Loader Service

(c) 2013 Schrodinger Inc.
'''

import Tkinter
from pymol import cmd, plugins, CmdException

def __init_plugin__(self=None):
    plugins.addmenuitem('PDB Loader Service', fetchdialog)

def get_trunc(var):
    return var.get().split(None, 1)[0]

def fetchdialog():
    from Tkinter import LEFT, RIGHT, TOP, BOTTOM

    app = plugins.get_pmgapp()
    root = plugins.get_tk_root()

    self = Tkinter.Toplevel(root)
    self.title('PDB Loader Service')
    self.minsize(250, 0)
    self.resizable(0,0)
    outer = self

    pad = 4

    type_options = [
        "pdb (Asymmetric Unit)",
        "pdb1 (Biological Unit)",
        "2fofc (Density)",
        "fofc (Difference Density)",
        "cid (PubChem Compound)",
        "sid (PubChem Substance)",
        "emd (EMDB Density)",
    ]

    var_code = Tkinter.StringVar(self)
    var_chain = Tkinter.StringVar(self)
    var_type = Tkinter.StringVar(self, type_options[0])
    var_keep = Tkinter.BooleanVar(self,
            not cmd.get_setting_boolean('autoclose_dialogs'))

    def callback(*args):
        code = var_code.get()
        type = get_trunc(var_type)
        if type == 'pdb':
            code += var_chain.get()
        try:
            result = cmd.fetch(code, type=type)
            if result == -1:
                raise CmdException('You entered an invalid pdb code: ' + code)
        except CmdException as e:
            import tkMessageBox
            tkMessageBox.showerror('Error', str(e), parent=self)
            return
        cmd.log('fetch %s, type=%s, async=0\n' % (code, type))
        if not var_keep.get():
            self.destroy()

    def callback_type(*args):
        v = get_trunc(var_type)
        if v.startswith('pdb') or v.endswith('fofc'):
            but_code.configure(width=4)
        else:
            but_code.configure(width=20)
        if v == 'pdb':
            frame_pdb.pack(side=LEFT)
        else:
            frame_pdb.pack_forget()

    def makerow(label=None, parent=None, **kwargs):
        master = Tkinter.Frame(parent or outer)
        master.pack(fill=Tkinter.X, **kwargs)
        if label:
            Tkinter.Label(master, text=label).pack(side=LEFT)
        return master

    padkw = {'padx': pad, 'pady': (pad, 0)}

    master = makerow("Type:", **padkw)
    but_type = Tkinter.OptionMenu(master, var_type, *type_options,
            command=callback_type).pack(side=LEFT)

    master = makerow("Code:", **padkw)
    but_code = Tkinter.Entry(master, textvariable=var_code, width=4)
    but_code.bind("<Return>", callback)
    but_code.pack(side=LEFT)

    frame_pdb = makerow("Chain:", master, side=LEFT, padx=8)
    but_chain = Tkinter.Entry(frame_pdb, textvariable=var_chain, width=4)
    but_chain.bind("<Return>", callback)
    but_chain.pack(side=LEFT)

    master = makerow(padx=pad, pady=(2*pad,0))
    but_ok = Tkinter.Button(master, width=10, text="OK",
            command=callback)
    but_cancel = Tkinter.Button(master, width=10, text="Cancel",
            command=self.destroy)

    but_cancel.pack(side=RIGHT, fill=Tkinter.X)
    but_ok.pack(side=RIGHT, fill=Tkinter.X)

    master = makerow(padx=pad, pady=(2,pad))
    but_keep = Tkinter.Checkbutton(master, text="Keep dialog open", variable=var_keep)
    but_keep.pack(side=RIGHT)

    but_code.focus_set()

# vi:expandtab:smarttab
