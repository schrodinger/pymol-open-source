
# this is quick and dirty irst crack at a gamess interface, written by someone who
# knows very little about the program (WLD) - hence the version 1 identifier...

# by the way, most of the below is untested...
      
import os
import shutil
import glob
import re
import string
import sys
import time

from chempy import feedback
from chempy.brick import Brick

atNum = {
   'H'  : 1,
   'C'  : 6,
   'N'  : 7,
   'O'  : 8,
   'F'  : 9,
   'P'  : 15,
   'S'  : 16,
   'Cl' : 17,
   'Br' : 35,
   'I'  : 53,
   }

# "do" is the preferred command for running tinker

def do(input,run_prefix='gamess_run',echo=None,
       punch=None,output=None):
   if feedback['gamess']:
      print " "+str(__name__)+': creating temporary files "%s.*"' % (run_prefix)
      print " "+str(__name__)+': launching gamess...' 
   try:
      for a in glob.glob(run_prefix+".*"):
         os.unlink(a)
   except:
      pass
   f = open(run_prefix+".inp",'w')
   for a in input:
      f.write(a)
   f.close()
   if echo:
      os.system(rungms_path+' '+run_prefix+" 2>&1 | tee "+run_prefix+".out")
   else:
      os.system(rungms_path+' '+run_prefix+" > "+run_prefix+".out 2>&1")      
# NFS workaround (flushes the directory cache so that glob will work)
   try: os.unlink(".sync")
   except: pass
   f = open(".sync",'w')
   f.close()
#
   if feedback['gamess']:
      print " "+str(__name__)+': job complete. '
   if punch:
      for src in glob.glob(run_prefix+".dat"):
         f = open(src)
         punch = f.readlines()
         f.close()
   if output:
      for src in glob.glob(run_prefix+".out"):
         f = open(src)
         output = f.readlines()
         f.close()
   return (output,punch)

if os.environ.has_key('GAMESS'):
   base = os.environ['GAMESS']
   bin_path = base + '/bin/'
   rungms_path = bin_path + 'rungms'
else:
   base = ''
   bin_path = ''
   params_path = ''

class State:

   def __init__(self):
      self.prefix = 'gamess_run'
      self.model = None
      self.data = None
      self.vec = None
      
   def load_model(self,model):
      self.model = model

   def get_data_group(self,basis = None,zmat = 1):
      model = self.model

      gmsList = []
      
      # write header records
      gmsList.append(" $DATA\n")
      gmsList.append(model.molecule.title+" from "+str(__name__)+"\n")
      gmsList.append("C1\n")

      # write atom records in an ordering compatible with internal
      # coordinate generation
      c = 1
      for z in model.get_internal_tuples():
         a = model.atom[z[0]]
         if not len(a.name):
            name = a.symbol + "%02d"%c
         else:
            name = a.name
         gmsList.append("%4s %5.2f %12.6f %12.6f %12.6f\n" %
                        (name,atNum[a.symbol],a.coord[0],
                         a.coord[1],a.coord[2]))
         c = c + 1
      gmsList.append(" $END\n")

      gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE. $END\n")
      return gmsList

   def get_contrl_group(self,runtyp='OPTIMIZE',exetyp='RUN',
                       coord='UNIQUE',nzvar = -1):
      gmsList = []
      model = self.model
      if nzvar:
         if nzvar<0:
            nzvar = (self.model.nAtom*3)-6
      gmsList.append(" $CONTRL RUNTYP=%s EXETYP=%s\n"
                     % (runtyp,exetyp) )
      if coord:
         gmsList.append("COORD=%s\n"%coord)
      if nzvar:
         gmsList.append("NZVAR=%d\n"%nzvar)
      chg = 0
      for a in model.atom:
         chg = chg + a.formal_charge
      chg = int(chg)
      if chg==0:
         icharg=None
      else:
         icharg='ICHARG=%d' % chg
      if icharg:
         gmsList.append("%s\n" % (icharg))
      gmsList.append(" $END\n")
      return gmsList

   def read_output_list(self,list):
      ll = len(list)
      c = 0
      crd_list = []
      chg_list = []
      for a in list:
         if a[0:36] == ' COORDINATES OF ALL ATOMS ARE (ANGS)':
            crd_list.append(c+3)
         if a[0:13] == ' NET CHARGES:':
            chg_list.append(c+4)
         c = c + 1
      atom = self.model.atom
      idx = {}
      c = 0
      for a in atom:
         idx[a.name]=c
         c = c + 1
      if len(crd_list):
         a = crd_list.pop()
         cc = 0
         while a<ll:
            l = list[a]
            name = string.strip(l[1:11])
            if name=='':
               break
            atom[idx[name]].coord = [float(l[16:31]),
                                     float(l[31:46]),
                                     float(l[46:61])]
            cc = cc + 1
            a = a + 1
         if cc and feedback['gamess']:
            print " "+str(__name__)+': coordinates modified for %d atoms.' % (cc)
      if len(chg_list):
         a = chg_list.pop()
         cc = 0
         while a<ll:
            l = list[a]
            name = string.strip(l[1:11])
            if name[0]=='-':
               break
            atom[idx[name]].partial_charge = float(l[19:27])
            a = a + 1
            cc = cc + 1
         if cc and feedback['gamess']:
            print " "+str(__name__)+': charges modified for %d atoms.' % (cc)
      
   def read_punch_list(self,list):
      ll = len(list)
      c = 0
      data_list = []
      vec_list = []
      for a in list:
         if a[0:6] == ' $DATA':
            data_list.append(c)
         elif a[0:5] == ' $VEC':
            vec_list.append(c)
         c = c + 1
      if len(data_list):
         a = data_list.pop()
         self.data = []
         data = self.data
         while a<ll:
            la = list[a]
            data.append(la)
            if la[0:5] == ' $END':
               break
            a = a + 1
         if feedback['gamess']:
            print " "+str(__name__)+': read new $DATA group.'
      if len(vec_list):
         a = vec_list.pop()
         self.vec = []
         vec = self.data
         while a<ll:
            la = list[a]
            vec.append(la)
            if la[0:5] == ' $END':
               break
            a = a + 1
         if feedback['gamess']:
            print " "+str(__name__)+': read new $VEC group.'

   def read_density_list(self,list,brick,z_step):
      ll = len(list)
      c = 0
      den_list = []
      for a in list:
         if a[0:37] == ' ELECTRON DENSITY, IPOINT,X,Y,Z,EDENS':
            den_list.append(c+1)
         c = c + 1
      if len(den_list):
         lst = 0
         a = den_list.pop()
         for x in xrange(brick.dim[0]):
            for y in xrange(brick.dim[1]):
               brick.lvl[x][y][z_step] = float(list[a][36:51])
               a = a + 1 
         if feedback['gamess']:
            print " "+str(__name__)+': read density slice %d of %d.' %(
               z_step+1,brick.dim[2])
         
   def get_basis_group(self,gbasis='N31',ngauss=6,ndfunc=1):
      gmsList = []
      gmsList.append(" $BASIS GBASIS=%s NGAUSS=%d NDFUNC=%d\n" %
                     (gbasis,ngauss,ndfunc))
      model = self.model
      chg = 0
      for a in model.atom:
         chg = chg + a.formal_charge
      chg = int(chg)
      if chg<0:
         diffsp=' DIFFSP=.TRUE.'
      else:
         diffsp=None
      if diffsp:
         gmsList.append("%s\n" % (diffsp))
      gmsList.append(" $END\n")
      return gmsList

   def get_zmat_group(self,auto=1,dlc=1):
      gmsList = []
      if auto and dlc:
         gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE. $END\n")
      else:
         raise RuntimeError
      return gmsList
   
   def get_eldens_group(self,morb=0):
      gmsList = []
      gmsList.append(" $ELDENS IEDEN=1 MORB=%i \n" % morb)
      gmsList.append("WHERE=GRID OUTPUT=PUNCH\n $END\n")
      return gmsList

   def get_grid_group(self,brick,z_step):

      origin = (
         brick.origin[0],
         brick.origin[1],
         brick.origin[2]+brick.grid[2]*z_step)
      x_coord = (
         brick.origin[0]+brick.range[0]+brick.grid[0]/100.0,
         brick.origin[1],
         brick.origin[2]+brick.grid[2]*z_step)
      y_coord = (
         brick.origin[0],
         brick.origin[1]+brick.range[1]+brick.grid[1]/100.0,
         brick.origin[2]+brick.grid[2]*z_step)
      gmsList = [
         " $GRID ORIGIN(1)=%12.5f,%12.5f,%12.5f\n" % origin,
         "XVEC(1) = %12.5f,%12.5f,%12.5f\n" % x_coord,
         "YVEC(1) = %12.5f,%12.5f,%12.5f\n" % y_coord,
         "SIZE = %12.5f\n" % brick.grid[0],
         " $END\n"
         ]
      return gmsList

   def get_scf(dirscf=1):
      gmsList = []
      if dirscf:
         gmsList.append(" $SCF DIRSCF=.TRUE. $END\n")
      return gmsList
   
   def get_optimize_job(self,dirscf=1):
      gmsList = []
      gmsList.extend(self.get_contrl_group())
      gmsList.extend(self.get_basis_group())
      gmsList.extend(self.get_scf())
      gmsList.extend(self.get_data_group())
      gmsList.extend(self.get_zmat_group())
      gmsList.append(" $STATPT NSTEP=50 $END\n")
      return gmsList
         
   def get_optimize_charge_job(self):
      gmsList = self.get_optimize_job()
      gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
      gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
      return gmsList

   def get_energy_charge_job(self):
      gmsList = self.get_energy_job()
      gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
      gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
      return gmsList

   def get_energy_job(self):
      gmsList=[]
      gmsList.extend(self.get_contrl_group(
         runtyp = 'ENERGY'
         ))
      gmsList.extend(self.get_basis_group())
      gmsList.extend(self.get_scf())
      gmsList.extend(self.get_data_group())
      gmsList.extend(self.get_zmat_group())
      return gmsList
   
   def get_density_job(self,brick,z_step):
      gmsList = []
      gmsList.extend(self.get_contrl_group(
         runtyp = 'PROP', 
         coord = None,
         nzvar = 0
         ))
      gmsList.extend(self.get_scf())
      gmsList.extend(self.data)
      gmsList.extend(self.vec)
      gmsList.extend(self.get_eldens_group())
      gmsList.extend(self.get_grid_group(brick,z_step))
      return gmsList
   
   def get_density(self,brick):
      for a in xrange(brick.dim[2]):
         gmsList = self.get_density_job(brick,a)
         result = do(gmsList,punch=1)
         self.read_density_list(result[1],brick,a)

   def get_charges(self):
      gmsList = self.get_energy_job()
      result = do(gmsList,output=1,punch=1)
      self.read_punch_list(result[1])
      
   def get_optimized_charges(self):
      gmsList = self.get_optimize_charge_job()
      result = do(gmsList,output=1,punch=1)
      self.read_density_list(result[1],brick,a)      
                   
if os.environ.has_key('GAMESS'):
   base = os.environ['GAMESS']
   rungms_path = base + '/rungms'
else:
   base = ''
   rungms_path = ''






