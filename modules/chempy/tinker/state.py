#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Scott Dixon, Metaphorics, LLC
#-* 
#-*
#Z* -------------------------------------------------------------------

from chempy import tinker,io,feedback
from chempy.tinker import keyword
import copy
import string

class State:

    def __init__(self):
        if feedback['verbose']:
            print ' '+str(self.__class__)+': created.'
        self.prefix = 'tinker'
        self.params = tinker.params_path+'chempy.prm'

        self.default = {}
        self.default['timestep'] = 1.0
        self.default['temperature'] = 300.0

        self.keywords = {}
        self.keywords['verbose'] = ''
        self.keywords['cutoff'] = 8.0      
        self.keywords['randomseed'] = 1234567890

        self.mapping = None
        self.echo = 0
        self.model = None
        self.reset_fragile()
        self.counter = 0

    def reset_fragile(self):
        self.restart = None
        self.frames = None
        
    def __write_keys(self,list):
        f=open(self.prefix+"_inp.key",'w')
        f.writelines(list)
        f.close()

    def analyze(self,kw=None,summary=1):
        if feedback['actions']:
            print ' '+str(self.__class__)+': starting energy run...'
        io.xyz.toFile(self.model,self.prefix+"_inp.xyz",
                          mapping = self.mapping)
        kw_list = [ "parameters "+self.params+"\n" ]
        for a in self.keywords.keys():
            kw_list.append("%s %s\n" % ( a,str(self.keywords[a])))
        kw_list.extend(keyword.get_partial_charge(self.model))
        if kw:
            kw_list.extend(kw)
        self.__write_keys(kw_list)
        tinker.do('analyze',
                     self.prefix+"_inp",
                     self.prefix+"_run",
                     self.prefix+"_out",
                     [ tinker.prefix,
                        'E',],
                     capture = 1 + self.echo )
        if summary:
# get energy information from output file
            f = open(self.prefix+"_out.out")
            self.summary = []
            flag = 0
            while 1:
                lin = f.readline()
                if not lin: break
                if not flag:
                    if lin[0:25]==' Total Potential Energy :':
                        self.summary.append([
                            string.strip(lin[0:23]),
                            float(lin[25:49])])
                        flag = 1
                else:
                    tok = string.split(string.strip(lin))
                    if len(tok):
                        if(tok[0]!='Energy'):
                            self.summary.append([
                                string.strip(lin[0:23]),
                                float(lin[25:49]),
                                int(lin[49:64])])
            f.close()
        else:
            self.summary = None
        if len(self.summary):
             self.energy = self.summary[0][1]
        else:
             self.energy = None

    def minimize(self,gradient=0.1,max_iter=100,
                     kw=None,summary=1):
        if feedback['actions']:
            print ' '+str(self.__class__)+': starting minimization run...'
        io.xyz.toFile(self.model,self.prefix+"_inp.xyz")
        self.counter = self.counter + max_iter
        kw_list = [ "parameters "+self.params+"\n" ]
        for a in self.keywords.keys():
            kw_list.append("%s %s\n" % ( a,str(self.keywords[a])))
        kw_list.extend(keyword.get_partial_charge(self.model))
        if kw:
            kw_list.extend(kw)
        kw_list.append("overwrite\n")
        kw_list.append("maxiter %d\n" % max_iter)
        self.__write_keys(kw_list)
        tinker.do('minimize',
                     self.prefix+"_inp",
                     self.prefix+"_run",
                     self.prefix+"_out",
                     [ tinker.prefix,
                        str(gradient),],
                     capture = 1 + self.echo )
        io.xyz.updateFromFile(self.model,self.prefix+"_out.xyz")
        self.restart = None
        if summary:
# get energy information from output file
            f = open(self.prefix+"_out.out")
            step = 1
            flag = 0
            self.summary = []
            while 1:
                lin = f.readline()
                if not lin: break
                if not flag:
                    if (lin[0:9]==' CG Iter ' or
                         lin[0:9]==' QN Iter '):
                        f.readline() # skip blank
                        flag = 1
                else:
                    try:
                        if int(lin[0:6])==step:
                            step = step + 1
                            self.summary.append([
                                int(lin[0:6]),
                                float(lin[6:19]),
                                float(lin[19:30]),
                                float(lin[30:41]),
                                float(lin[41:51]),
                                float(lin[51:60]),
                                int(lin[60:67]),
                                string.strip(lin[67:])])
                    except ValueError:
                        pass
            f.close()
        else:
            self.summary = None
        if len(self.summary):
             self.energy = self.summary[-1][1]
        else:
             self.energy = None
            
    def dynamics(self,steps=100,
                     timestep=None,
                     interval=None,
                     temperature=None,
                     kw=None,summary=1):
        if feedback['actions']:
            print ' '+str(self.__class__)+': starting dynamics run...'
        io.xyz.toFile(self.model,self.prefix+"_inp.xyz")
        self.counter = self.counter + steps * timestep
# if restart is available, then use it.
        if self.restart:
            f = open(self.prefix+"_inp.dyn","w")
            f.writelines(self.restart)
            f.close()
            if feedback['actions']:
                print ' '+str(self.__class__)+': using restart information...'
#         
        kw_list = [ "parameters "+self.params+"\n" ]
        for a in self.keywords.keys():
            kw_list.append("%s %s\n" % ( a,str(self.keywords[a])))
        kw_list.extend(keyword.get_partial_charge(self.model))
        if kw:
            kw_list.extend(kw)
        kw_list.append("archive\n")
        self.__write_keys(kw_list)
        if not timestep:
            timestep = self.default['timestep']
        if not temperature:
            temperature = self.default['temperature']
        if not interval:
            interval = (timestep*steps)/1000.0
        else:
            interval = interval/1000.0
        tinker.do('dynamic',
                     self.prefix+"_inp",
                     self.prefix+"_run",
                     self.prefix+"_out",
                     [ tinker.prefix,
                        str(steps),
                        str(timestep),
                        str(interval),
                        str(temperature)],
                     capture = 1 + self.echo )
        self.frames = io.arc.fromFile(self.prefix+"_out.arc")
        if len(self.frames):
            c = 0
            atm = self.model.atom
            for a in self.frames[-1]:
                atm[c].coord = copy.deepcopy(a)
                c = c + 1
# get restart information
        f = open(self.prefix+"_out.dyn")
        self.restart = f.readlines()
        f.close()
# get energy information from output file
        if summary:
            f = open(self.prefix+"_out.out")
            step = 1
            flag = 0
            self.summary = []
            while 1:
                lin = f.readline()
                if not lin: break
                if not flag:
                    if lin[0:12]=='    MD Step ':
                        f.readline() # skip blank
                        flag = 1
                else:
                    try:
                        if int(lin[0:10])==step:
                            step = step + 1
                            self.summary.append([
                                int(lin[0:10]),
                                float(lin[10:24]),
                                float(lin[24:38]),
                                float(lin[38:50]),
                                float(lin[50:61])])
                    except ValueError:
                        pass
            f.close()
        else:
            self.summary = None
        if len(self.summary):
            self.energy = self.summary[-1][1]
        else:
            self.energy = None

        
    def load_model(self,a):
        if feedback['verbose']:
            print ' '+str(self.__class__)+': new model loaded.'
        self.model = a
        self.reset_fragile()





