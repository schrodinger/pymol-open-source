'''
Simple Text Editor
'''

import os
import sys
import Tkinter
import tkFileDialog
import tkMessageBox

class TextEditor:

    def _write(self, handle):
        content = self._get()
        handle.write(content)
        self._savedcontent = content

    def _open(self, filename):
        self.filename = filename or ''
        if filename and os.path.exists(filename):
            with open(filename, 'rU') as handle:
                content = handle.read()
        else:
            content = ''
        self._set(content)

    def _get(self):
        return self.text.get(0.0, Tkinter.END)

    def _set(self, content):
        self.text.delete(0.0, Tkinter.END)
        self.text.insert(0.0, content)
        self._savedcontent = self._get()

    def doSaveAs(self, *args):
        handle = tkFileDialog.asksaveasfile(mode='w',
                initialfile=os.path.basename(self.filename),
                initialdir=os.path.dirname(self.filename),
                parent=self.root)
        if handle:
            with handle:
                self._write(handle)
                self.filename = handle.name

    def doSave(self, *args):
        if not self.filename:
            return self.doSaveAs()
        with open(self.filename, 'w') as handle:
            self._write(handle)

    def doOpen(self, *args):
        filename = tkFileDialog.askopenfilename(parent=self.root)
        if filename:
            self._open(filename)

    def onClose(self):
        if self._get() != self._savedcontent:
            ok = tkMessageBox.askyesnocancel("Save?", "Save before quit?",
                    parent=self.root)
            if ok:
                self.doSave()
            elif ok is None:
                return
        self.root.destroy()

    def __init__(self, parent=None, filename='', title='Text Editor'):
        self.root = Tkinter.Toplevel(parent) if parent else Tkinter.Tk()
        self.root.title(title)
        self.root.minsize(width=500, height=400)
        self.root.protocol("WM_DELETE_WINDOW", self.onClose)

        menubar = Tkinter.Menu(self.root)
        filemenu = Tkinter.Menu(menubar)
        filemenu.add_command(label="Open", command=self.doOpen, accelerator="Ctrl+O")
        filemenu.add_command(label="Save", command=self.doSave, accelerator="Ctrl+S")
        filemenu.add_command(label="Save as ...", command=self.doSaveAs, accelerator="Ctrl+Shift+S")
        menubar.add_cascade(label="File", menu=filemenu)
        self.root.config(menu=menubar)

        self.text = Tkinter.Text(self.root, background='white', foreground='black')
        self.text.pack(expand=Tkinter.YES, fill=Tkinter.BOTH)

        self._open(filename)

        self.text.bind("<Control-o>", self.doOpen)
        self.text.bind("<Control-s>", self.doSave)
        self.text.bind("<Control-S>", self.doSaveAs)

if __name__ == '__main__':
    try:
        filename = sys.argv[1]
    except:
        filename = ''
    app = TextEditor(None, filename)
    app.root.mainloop()
