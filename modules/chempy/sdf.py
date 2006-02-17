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

class SDFRec:
    
    def __init__(self,sdflist):
        getkee = re.compile("^>\s+<([^>]*)>")
        gettag = re.compile("^>\s+<[^>]*>\s+\((.*)\)")
        ll = len(sdflist)
        if ll<4:
            print " SDFRec: invalid SDF record format #1"
            raise RuntimeError
        self.kees = ['MOL'] # separate key list to preserve order
        self.data = {}
        self.ref_code = {}
        self.data['MOL'] = []
        mol = self.data['MOL']
        l = 0
        while l<ll:
            mol.append(sdflist[l])
            if (sdflist[l][0:6]=='M  END') or (sdflist[l][0:5]=='M END'):
                break;
            if sdflist[l][0:1]=='>':
                mol[len(mol)-1]='M  END\n'
                sdflist.insert(l,'M  END\n')
                ll=len(sdflist)
                break;
            l = l + 1
            if l>=ll:
                print " SDFRec: invalid SDF record format #2"
                raise RuntimeError
        while l<ll:
            if sdflist[l][0]=='>':
                sl = sdflist[l]
                kee_match = getkee.match(sl)
                if not kee_match:
                    print " SDFRec: invalid SDF record format #3"
                    raise RuntimeError
                kee = kee_match.group(1)
                self.kees.append(kee)
                ref_code_match = gettag.match(sl)
                if ref_code_match:
                    self.ref_code[kee] = ref_code_match.group(1)
                else:
                    self.ref_code[kee] = ''
                self.data[kee] = []
                sd = self.data[kee]
                l = l + 1
                while l<ll:
                    if len(string.strip(sdflist[l]))!=0:
                        sd.append(sdflist[l])
                        l = l + 1
                    else:
                        break;
            else:
                l = l + 1

    def toList(self):
        r = []
        for k in self.kees:
            if k!='MOL':
                if self.ref_code[k]!='':
                    r.append(">  <"+k+"> ("+self.ref_code[k]+")\n")
                else:
                    r.append(">  <"+k+">\n")               
            for a in self.data[k]:
                r.append(a)
            if k!='MOL':
                r.append("\n")
        return r

    def get(self,kee):
        if self.data.has_key(kee):
            return self.data[kee]
        else:
            return None

    def get_single(self,kee): # automatic stripping
        if self.data.has_key(kee):
            sdk = self.data[kee]
            if len(sdk):
                return string.strip(sdk[0])
            else:
                return None
        else:
            return None
            
    def set_single(self,kee,data,ref_code=None): # adds LF
        self.set(kee,[data+'\n'],ref_code)
        
    def get_model(self):
        return io.mol.fromList(self.get('MOL'))
        
    def set_model(self,model):
        self.set('MOL',io.mol.toList(model))
                    
    def set(self,kee,data,ref_code=None):
        if kee not in self.kees:
            self.kees.append(kee)
            self.ref_code[kee]=''
        if ref_code!=None:
            self.ref_code[kee]=ref_code
        self.data[kee] = copy.deepcopy(data)

    def delete(self,kee):
        self.kees.remove(kee)
        del self.data[kee]

class SDF:
    
    def __init__(*args):
        mode = 'r'
        if len(args)<2:
            raise ValueError
        self = args[0]
        fname = args[1]
        if len(args)==3:
            mode = args[2]
        self.mode = mode
        self.at_eof = 0
        if mode not in ('w','r','wa','pf','url'):
            print " SDF: bad mode"
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
        
    def read(self): # returns SDFRec or None at end of file
        cur = []
        while 1:
            s = self.file.readline()
            if not s:
                return None
            elif s[0:4]==r'$$$$':
                return SDFRec(cur)
            else:
                cur.append(s)
            
    def close(self):
        self.file.close()
        
    
