
import string
import os


#
# Basic chempy types
#

class Atom:

   defaults = {
      'symbol'              : 'X',
      'name'                : '',
      'resn'                : 'UNK',
      'resn_code'           : 'X',
      'resi'                : '1',
      'resi_number'         : 1,
      'b'                   : 0.0,
      'q'                   : 1.0,
      'vdw'                 : 0.0,
      'alt'                 : '',
      'hetatm'              : 1,
      'segi'                : '',
      'chain'               : '',
      'coord'               : [9999.999,9999.999,9999.999],
      'formal_charge'       : 0.0,
      'partial_charge'      : 0.0,
# Flags
      'flags'               : 0,
# Force-fields
      'numeric_type'        : -9999,
      'text_type'           : '??',
# MDL Mol-files
      'stereo'              : 0,
# Macromodel files
      'color_code'          : 2,
      }
   
   def __getattr__(self,attr):
      if Atom.defaults.has_key(attr):
         return Atom.defaults[attr]
      else:
         raise AttributeError(attr)

   def get_mass(self):
      return atomic_mass[self.symbol]
   
   def has(self,attr):
      return self.__dict__.has_key(attr) 

   def in_same_residue(self,other):
      if self.resi == other.resi:
         if self.chain == other.chain:
            if self.segi == other.segi:
               return 1
      return 0

   def new_in_residue(self):
      newat = Atom()
      if self.has('segi'):        newat.segi        = self.segi
      if self.has('chain'):       newat.chain       = self.chain
      if self.has('resn'):        newat.resn        = self.resn
      if self.has('resn_code'):   newat.resn_code   = self.resn_code
      if self.has('resi'):        newat.resi        = self.resi
      if self.has('resi_number'): newat.resi_number = self.resi_number
      if self.has('hetatm'):      newat.hetatm      = self.hetatm
      return newat

   def get_signature(self):
      return string.join([self.segi,self.chain,self.resn,
                          self.resi,self.symbol,self.name],':')
   
   def __cmp__(self,other):
      if type(self)==type(other):
         if self.segi == other.segi:
            if self.chain == other.chain:
               if self.resi_number == other.resi_number:
                  if self.resn == other.resn:
                     if self.resi == other.resi:
                        if self.symbol == other.symbol:
                           if self.name == other.name:
                              return cmp(id(self),id(other))
                           else:
                              return cmp(self.name,other.name)
                        else:
                           return cmp(self.symbol,other.symbol)
                     else:
                        return cmp(self.resi,other.resi)
                  else:
                     return cmp(self.resn,other.resn)
               else:
                  return cmp(self.resi_number,other.resi_number)               
            else:
               return cmp(self.chain,other.chain)
         else:
            return cmp(self.segi,other.segi)
      else:
         return cmp(type(self),type(other))
      
class Bond:

   defaults = {
      'order'           : 1,
      'stereo'          : 0
      }

   def __getattr__(self,attr):
      if Bond.defaults.has_key(attr):
         return Bond.defaults[attr]
      else:
         raise AttributeError(attr)
      
   def has(self,attr):
      return self.__dict__.has_key(attr) 

class Molecule:

   defaults = {
      'dim_code'        : '3D',
      'title'           : 'untitled',
      'comments'        : '',
      'chiral'          : 1
      }

   def __getattr__(self,attr):
      if Molecule.defaults.has_key(attr):
         return Molecule.defaults[attr]
      else:
         raise AttributeError(attr)
      
   def has(self,attr):
      return self.__dict__.has_key(attr) 
   
class Storage:

   def updateFromList(self,indexed,**params):
      pass
   
   def fromList(self,**params):
      return chempy.indexed()
   
   def toList(self,indexed,**params):
      return []

   def updateFromFile(self,indexed,fname,**params):
      fp = open(fname)
      result = apply(self.updateFromList,(indexed,fp.readlines()),params)
      fp.close()

   def fromFile(self,fname,**params):
      if feedback['io']:
         print ' chempy: reading "%s".' % fname
      fp = open(fname)
      result = apply(self.fromList,(fp.readlines(),),params)
      fp.close()
      return result

   def toFile(self,indexed,fname,**params):
      if feedback['io']:
         print ' chempy: writing "%s".' % fname
      fp = open(fname,'w')
      result = fp.writelines(apply(self.toList,(indexed,),params))
      fp.close()

feedback = { 'warnings': 1,
             'terse'   : 1,
             'io'      : 1,
             'actions' : 1,
             'tinker'  : 1,
             'gamess'  : 1,             
             'atoms'   : 0,
             'bonds'   : 0,                          
             'verbose' : 0,
             }

if os.environ.has_key('CHEMPY_PATH'):
   path = os.environ['CHEMPY_PATH'] + '/'
elif os.environ.has_key('FREEMOL_MODULES'):
   path = os.environ['FREEMOL_MODULES'] + '/chempy/'
elif os.environ.has_key('PYMOL_PATH'):
   path = os.environ['PYMOL_PATH'] + '/modules/chempy/'
else:
   path = ''

# double check these values...

atomic_mass = {
   'H'  :   1.008,
   'C'  :  12.011,
   'N'  :  14.006,
   'O'  :  15.999,
   'F'  :  18.998,
   'Cl' :  35.453,
   'Br' :  79.904,
   'I'  : 126.904,
   'S'  :  32.064,
   'Na' :  22.990,
   'K'  :  39.102,
   'Cu' :  63.546,
   'Zn' :  65.370,
   'Mg' :  24.312,
   'Ca' :  40.080,
   'P'  :  30.974
   }

