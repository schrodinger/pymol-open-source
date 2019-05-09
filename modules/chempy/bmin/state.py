from __future__ import print_function

from chempy import bmin,feedback
from chempy import io

import os
import re
import getpass # for getuser()

class State:

    def __init__(self):
        if feedback['verbose']:
            print(' '+str(self.__class__)+': created.')
        self.default = {}
        self.echo = 0
        self.model = None
        self.counter = 0
        self.prefix = "bmintmp"

    def minimize(self,max_iter=100,fix_flag=None,rest_flag=None,
                     rest_coeff = 100.0,solvation=None):
        if feedback['actions']:
            print(' '+str(self.__class__)+': starting minimization run...')
        io.mmd.toFile(self.model,self.prefix+".dat")

        f = open(self.prefix+".com",'w')
        # get home-relative path
#      pth = os.getcwd()
#      pth = re.sub(r".*\/"+getpass.getuser()+"\/",'',pth)
        # provide io filenames
#      f.write("%s\n%s\n"%(pth+"/"+self.prefix+".dat",
#                          pth+"/"+self.prefix+".out"))
        f.write("%s\n%s\n"%(self.prefix+".dat",
                                  self.prefix+".out"))
        f.write(" MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n")
        # select forcefield treatments
        if not solvation: # no solvent, constant dielectric
            f.write(" FFLD      10      1      0      1     1.0000     0.0000     0.0000     0.0000\n")
        else:
            f.write(''' FFLD      10      1      0      1     1.0000     0.0000     0.0000     0.0000
 SOLV       3      1      0      0     0.0000     0.0000     0.0000     0.0000
 EXNB       0      0      0      0     0.0000     0.0000     0.0000     0.0000
''')
        # read files
        f.write(" READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n")
        # fix/restrain atoms according to flags provided
        if fix_flag is not None: # are we fixing any atoms?
            c = 0
            mask= 2 ** fix_flag
            for a in self.model.atom:
                c = c + 1
                if mask&a.flags:
                    f.write(" FXAT  %6d      0      0      0    -1.0000     0.0000     0.0000     0.0000\n"%
                              c)
        if rest_flag is not None: # are we restraining any atoms?
            c = 0
            mask= 2 ** rest_flag
            for a in self.model.atom:
                c = c + 1
                if mask&a.flags:
                    f.write(" FXAT  %6d      0      0      0 %10.4f     0.0000     0.0000     0.0000\n"%
                              (c,rest_coeff))
        f.write(
''' CONV       2      0      0      0     0.0500     0.0000     0.0000     0.0000
 MINI       1      0 %6d      0     0.0000     0.0000     0.0000     0.0000
 DEBG 6
 '''%(max_iter))
        f.close()

        bmin.do(self.prefix)
        io.mmd.updateFromFile(self.model,self.prefix+".out")
        if hasattr(self.model.molecule,'energy'):
            self.model.molecule.title = "%1.3f"%self.model.molecule.energy
    def load_model(self,a):
        if feedback['verbose']:
            print(' '+str(self.__class__)+': new model loaded.')
        self.model = a
