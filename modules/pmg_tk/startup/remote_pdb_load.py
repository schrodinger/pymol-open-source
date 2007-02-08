# Copyright Notice
# ================
# 
# The PyMOL Plugin source code in this file is copyrighted, but you can
# freely use and copy it as long as you don't change or remove any of
# the copyright notices.
# 
# ----------------------------------------------------------------------
# This PyMOL Plugin is Copyright (C) 2004 by Charles Moad <cmoad@indiana.edu>
# 
#                        All Rights Reserved
# 
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name(s) of the author(s) not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
# 
# THE AUTHOR(S) DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
# INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS.  IN
# NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY SPECIAL, INDIRECT OR
# CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
# USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
# OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

from Tkinter import *
from pymol import cmd

def __init__(self):
    # Simply add the menu entry and callback
    self.menuBar.addmenuitem('Plugin', 'command',
                                     'PDB Loader Service',
                                     label = 'PDB Loader Service',
                                     command = lambda s=self : FetchPDB(s))
#    FetchPDB(self)
    
# Class that simply takes the pdbcode from a dialog and retrieves the file
class FetchPDB:
    def __init__(self, app):
        import tkSimpleDialog
        import tkMessageBox
        import urllib
        import gzip
        import os
        import string

        pdbCode = tkSimpleDialog.askstring('PDB Loader Service',
                                                      'Please enter a 4-digit pdb code:',
                                                      parent=app.root)
        
        if pdbCode: # None is returned for user cancel
            pdbCode = string.upper(pdbCode)
            try:
#                filename = urllib.urlretrieve('http://www.rcsb.org/pdb/cgi/export.cgi/' +
#                                                        pdbCode + '.pdb.gz?format=PDB&pdbId=' +
#                                                        pdbCode + '&compression=gz')[0]
                filename = urllib.urlretrieve('http://www.rcsb.org/pdb/files/'
                                              + pdbCode + '.pdb.gz')[0]
            except:
                tkMessageBox.showerror('Connection Error',
                                       'Can not access to the PDB database.\n'+
                                       'Please check your Internet access.',
                                       parent=app.root)
            else:
                if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
                    # Uncompress the file while reading
                    fpin = gzip.open(filename)

                    # Form the pdb output name
                    outputname = os.path.dirname(filename) + os.sep + pdbCode + '.pdb'
                    fpout = open(outputname, 'w')
                    fpout.write(fpin.read()) # Write pdb file

                    fpin.close()
                    fpout.close()

                    cmd.load(outputname,quiet=0) # Load the fresh pdb
                else:
                    tkMessageBox.showerror('Invalid Code',
                                                  'You entered an invalid pdb code:' + pdbCode,
                                                  parent=app.root)

                os.remove(filename) # Remove tmp file (leave the pdb)

