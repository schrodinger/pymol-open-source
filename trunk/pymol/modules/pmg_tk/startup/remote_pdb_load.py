#
# This simple plugin will automatically load a pdb from the
# web given the pdb code.  It does error handling as well.
#
# Charles Moad (cmoad@indiana.edu)
# Scientific Data Analysis Lab
# Indiana University
#

from Tkinter import *
from pymol import cmd

def __init__(self):
   # Simply add the menu entry and callback
   self.menuBar.addmenuitem('Plugin', 'command',
                            'PDB Loader Service',
                            label = 'PDB Loader Service',
                            command = lambda s=self : FetchPDB(s))

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
         pdcCode = string.upper(pdbCode)
         filename = urllib.urlretrieve('http://www.rcsb.org/pdb/cgi/export.cgi/' +
                                       pdbCode + '.pdb.gz?format=PDB&pdbId=' +
                                       pdbCode + '&compression=gz')[0]
         if (os.path.getsize(filename) > 0): # If 0, then pdb code was invalid
            # Uncompress the file while reading
            fpin = gzip.open(filename)

            # Form the pdb output name
            outputname = os.path.dirname(filename) + os.sep + pdbCode + '.pdb'
            fpout = open(outputname, 'w')
            fpout.write(fpin.read()) # Write pdb file
            
            fpin.close()
            fpout.close()
            
            cmd.load(outputname) # Load the fresh pdb

         else:
            tkMessageBox.showerror('Invalid Code',
                                   'You entered an invalid pdb code:' + pdbCode,
                                   parent=app.root)

         os.remove(filename) # Remove tmp file (leave the pdb)
