#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

import string
import re
import copy

from chempy import io

class CIFRec:
    
    def __init__(self,cif_list):
        print cif_list
    
    def toList(self):
        return r

class CIF:
    
    def __init__(*args):
        mode = 'r'
        if len(args)<2:
            raise ValueError
        self = args[0]
        self.input_line = None
        fname = args[1]
        if len(args)==3:
            mode = args[2]
        self.mode = mode
        self.at_eof = 0
        if mode not in ('w','r','wa','pf','url'):
            print " CIF: bad mode"
            return None
        if mode=='pf': # pseudofile
            self.file = fname
        elif (mode[0:1]=='r') and (string.find(fname,':')>1):
            # does this look like a URL? (but not a DOS path)
            from urllib import urlopen
            self.file = urlopen(fname)
        else:
            self.file = open(fname,mode)

    def write(self,rec):
        lst = rec.toList()
        for a in lst:
            self.file.write(a)
        self.file.write('$$$$\n')
        
    def read(self): # returns CIFRec or None at end of file
        cur = []
        while 1:
            if self.input_line == None:
                s = self.file.readline()
            else:
                s = self.input_line
                self.input_line = None
            if not s:
                return None
            elif s[0:5] == r'data_':
                if len(cur)>0:
                    self.input_line = s
                    return CIFRec(cur)
                else:
                    cur.append(s)
            else:
                cur.append(s)
            
    def close(self):
        self.file.close()
        
    
