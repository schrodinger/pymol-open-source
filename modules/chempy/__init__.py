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
import os
import copy

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
# Secondary structure
        'ss'                  : '',
        }
    
    def __getattr__(self,attr):
        if Atom.defaults.has_key(attr):
            return copy.deepcopy(Atom.defaults[attr])
        else:
            raise AttributeError(attr)

    def get_mass(self):
        '''Given the chemical symbol the atomic mass is returned'''      
        return atomic_mass[self.symbol]

    def get_number(self):
        '''Given the chemical symbol the atomic number is returned'''
        return atomic_number[self.symbol]

    def get_implicit_valence(self):
        return implicit_valence[self.symbol]
    
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
        'chiral'          : 1,
        'spacegroup'      : 'P 1',
        'cell'            : [1.0, 1.0, 1.0, 90.0, 90.0, 90.0],
        }

    def __getattr__(self,attr):
        if Molecule.defaults.has_key(attr):
            return Molecule.defaults[attr]
        else:
            raise AttributeError(attr)
        
    def has(self,attr):
        return self.__dict__.has_key(attr) 
    
class Storage:

    def my_open(self,fname,mode='r'):
        if (mode[0:1]=='r') and (string.find(fname,':')>1):
            import urllib
            return urllib.urlopen(fname)
        else:
            return open(fname,mode)
        
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
        fp = self.my_open(fname)
        result = apply(self.fromList,(fp.readlines(),),params)
        fp.close()
        return result

    def toFile(self,indexed,fname,**params):
        if feedback['io']:
            print ' chempy: writing "%s".' % fname
        fp = open(fname,'w')
        result = fp.writelines(apply(self.toList,(indexed,),params))
        fp.close()

class PseudoFile:

    def __init__(self,list=[]):
        self.list = copy.deepcopy(list)

    def write(self,st):
        self.list.append(str(st))
    
    def readline(self):
        try:
            return self.list.pop(0)
        except:
            return None

    def close(self):
        self.list = None
  
feedback = { 'warnings': 1,
                 'terse'   : 1,
                 'io'      : 1,
                 'actions' : 1,
                 'tinker'  : 1,
                 'gamess'  : 1,             
                 'atoms'   : 0,
                 'bonds'   : 0,                          
                 'verbose' : 0,
                 'bmin'    : 1,
                 }

if os.environ.has_key('CHEMPY_DATA'):  # 
    path = os.environ['CHEMPY_DATA'] + '/'
elif os.environ.has_key('PYMOL_DATA'):
    path = os.environ['PYMOL_DATA'] + '/chempy/'
elif os.environ.has_key('PYMOL_PATH'):
    path = os.environ['PYMOL_PATH'] + '/data/chempy/'   
elif os.environ.has_key('FREEMOL_MODULES'):
    path = os.environ['FREEMOL_MODULES'] + '/chempy/'
else:
    path = ''

# double check these values...
#hvd values obtained from http://www.webelements.com/ and recorded to their
#    known accuracy.

atomic_mass = {
    'H'  :   1.00794,
    'He' :   4.002602,
    'HE' :   4.002602,
    'Li' :   6.941,
    'LI' :   6.941,
    'Be' :   9.012182,
    'BE' :   9.012182,
    'B'  :  10.811,
    'C'  :  12.0107,
    'N'  :  14.0067,
    'O'  :  15.9994,
    'F'  :  18.9984032,
    'Ne' :  20.1797,
    'NE' :  20.1797,
    'Na' :  22.989770,
    'NA' :  22.989770,
    'Mg' :  24.3050,
    'MG' :  24.3050,
    'Al' :  26.981538,
    'AL' :  26.981538,
    'Si' :  28.0855,
    'SI' :  28.0855,
    'P'  :  30.973761,
    'S'  :  32.065,
    'Cl' :  35.453,
    'CL' :  35.453,
    'Ar' :  39.948,
    'AR' :  39.948,
    'K'  :  39.0983,
    'Ca' :  40.078,
    'CA' :  40.078,
    'Sc' :  44.955910,
    'SC' :  44.955910,
    'Ti' :  47.867,
    'TI' :  47.867,
    'V'  :  50.9415,
    'Cr' :  51.9961,
    'CR' :  51.9961,
    'Mn' :  54.938049,
    'MN' :  54.938049,
    'Fe' :  55.845,
    'FE' :  55.845,
    'Co' :  58.933200,
    'CO' :  58.933200,
    'Ni' :  58.6934,
    'NI' :  58.6934,
    'Cu' :  63.546,
    'CU' :  63.546,
    'Zn' :  65.39,
    'ZN' :  65.39,
    'Ga' :  69.723,
    'GA' :  69.723,
    'Ge' :  72.64,
    'GE' :  72.64,
    'As' :  74.92160,
    'AS' :  74.92160,
    'Se' :  78.96,
    'SE' :  78.96,
    'Br' :  79.904,
    'BR' :  79.904,   
    'Kr' :  83.80,
    'KR' :  83.80,
    'Rb' :  85.4678,
    'RB' :  85.4678,
    'Sr' :  87.62,
    'SR' :  87.62,
    'Y'  :  88.90585,
    'Zr' :  91.224,
    'ZR' :  91.224,
    'Nb' :  92.90638,
    'NB' :  92.90638,
    'Mo' :  95.94,
    'MO' :  95.94,
    'Tc' :  98,
    'TC' :  98,
    'Ru' : 101.07,
    'RU' : 101.07,
    'Rh' : 102.90550,
    'RH' : 102.90550,
    'Pd' : 106.42,
    'PD' : 106.42,
    'Ag' : 107.8682,
    'AG' : 107.8682,
    'Cd' : 112.411,
    'CD' : 112.411,
    'In' : 114.818,
    'IN' : 114.818,
    'Sn' : 118.710,
    'SN' : 118.710,
    'Sb' : 121.760,
    'SB' : 121.760,
    'Te' : 127.60,
    'TE' : 127.60,
    'I'  : 126.90447,
    'Xe' : 131.293,
    'XE' : 131.293,
    'Cs' : 132.90545,
    'CS' : 132.90545,
    'Ba' : 137.327,
    'BA' : 137.327,
    'La' : 138.9055,
    'LA' : 138.9055,
    'Ce' : 140.116,
    'CE' : 140.116,
    'Pr' : 140.90765,
    'PR' : 140.90765,
    'Nd' : 144.24,
    'ND' : 144.24,
    'Pm' : 145,
    'PM' : 145,
    'Sm' : 150.36,
    'SM' : 150.36,
    'Eu' : 151.964,
    'EU' : 151.964,
    'Gd' : 157.25,
    'GD' : 157.25,
    'Tb' : 158.92534,
    'TB' : 158.92534,
    'Dy' : 162.50,
    'DY' : 162.50,
    'Ho' : 164.93032,
    'HO' : 164.93032,
    'Er' : 167.259,
    'ER' : 167.259,
    'Tm' : 168.93421,
    'TM' : 168.93421,
    'Yb' : 173.04,
    'YB' : 173.04,
    'Lu' : 174.967,
    'LU' : 174.967,
    'Hf' : 178.49,
    'HF' : 178.49,
    'Ta' : 180.9479,
    'TA' : 180.9479,
    'W'  : 183.84,
    'Re' : 186.207,
    'RE' : 186.207,
    'Os' : 190.23,
    'OS' : 190.23,
    'Ir' : 192.217,
    'IR' : 192.217,
    'Pt' : 195.078,
    'PT' : 195.078,
    'Au' : 196.96655,
    'AU' : 196.96655,
    'Hg' : 200.59,
    'HG' : 200.59,
    'Tl' : 204.3833,
    'TL' : 204.3833,
    'Pb' : 207.2,
    'PB' : 207.2,
    'Bi' : 208.98038,
    'BI' : 208.98038,
    'Po' : 208.98,
    'PO' : 208.98,
    'At' : 209.99,
    'AT' : 209.99,
    'Rn' : 222.02,
    'RN' : 222.02,
    'Fr' : 223.02,
    'FR' : 223.02,
    'Ra' : 226.03,
    'RA' : 226.03,
    'Ac' : 227.03,
    'AC' : 227.03,
    'Th' : 232.0381,
    'TH' : 232.0381,
    'Pa' : 231.03588,
    'PA' : 231.03588,
    'U'  : 238.02891,
    'Np' : 237.05,
    'NP' : 237.05,
    'Pu' : 244.06,
    'PU' : 244.06,
    'Am' : 243.06,
    'AM' : 243.06,
    'Cm' : 247.07,
    'CM' : 247.07,
    'Bk' : 247.07,
    'BK' : 247.07,
    'Cf' : 251.08,
    'CF' : 251.08,
    'Es' : 252.08,
    'ES' : 252.08,
    'Fm' : 257.10,
    'FM' : 257.10,
    'Md' : 258.10,
    'MD' : 258.10,
    'No' : 259.10,
    'NO' : 259.10,
    'Lr' : 262.11,
    'LR' : 262.11,
    'Rf' : 261.11,
    'RF' : 261.11,
    'Db' : 262.11,
    'DB' : 262.11,
    'Sg' : 266.12,
    'SG' : 266.12,
    'Bh' : 264.12,
    'BH' : 264.12,
    'Hs' : 269.13,
    'HS' : 269.13,
    'Mt' : 268.14,
    'MT' : 268.14,
    }

atomic_number = {
    'H'  :   1,
    'He' :   2,
    'HE' :   2,
    'Li' :   3,
    'LI' :   3,
    'Be' :   4,
    'BE' :   4,
    'B'  :   5,
    'C'  :   6,
    'N'  :   7,
    'O'  :   8,
    'F'  :   9,
    'Ne' :  10,
    'NE' :  10,
    'Na' :  11,
    'NA' :  11,
    'Mg' :  12,
    'MG' :  12,
    'Al' :  13,
    'AL' :  13,
    'Si' :  14,
    'SI' :  14,
    'P'  :  15,
    'S'  :  16,
    'Cl' :  17,
    'CL' :  17,
    'Ar' :  18,
    'AR' :  18,
    'K'  :  19,
    'Ca' :  20,
    'CA' :  20,
    'Sc' :  21,
    'SC' :  21,
    'Ti' :  22,
    'TI' :  22,
    'V'  :  23,
    'Cr' :  24,
    'CR' :  24,
    'Mn' :  25,
    'MN' :  25,
    'Fe' :  26,
    'FE' :  26,
    'Co' :  27,
    'CO' :  27,
    'Ni' :  28,
    'NI' :  28,
    'Cu' :  29,
    'CU' :  29,
    'Zn' :  30,
    'ZN' :  30,
    'Ga' :  31,
    'GA' :  31,
    'Ge' :  32,
    'GE' :  32,
    'As' :  33,
    'AS' :  33,
    'Se' :  34,
    'SE' :  34,
    'Br' :  35,
    'BR' :  35,
    'Kr' :  36,
    'KR' :  36,
    'Rb' :  37,
    'RB' :  37,
    'Sr' :  38,
    'SR' :  38,
    'Y'  :  39,
    'Zr' :  40,
    'ZR' :  40,
    'Nb' :  41,
    'NB' :  41,
    'Mo' :  42,
    'MO' :  42,
    'Tc' :  43,
    'TC' :  43,
    'Ru' :  44,
    'RU' :  44,
    'Rh' :  45,
    'RH' :  45,
    'Pd' :  46,
    'PD' :  46,
    'Ag' :  47,
    'AG' :  47,
    'Cd' :  48,
    'CD' :  48,
    'In' :  49,
    'IN' :  49,
    'Sn' :  50,
    'SN' :  50,
    'Sb' :  51,
    'SB' :  51,
    'Te' :  52,
    'TE' :  52,
    'I'  :  53,
    'Xe' :  54,
    'XE' :  54,
    'Cs' :  55,
    'CS' :  55,
    'Ba' :  56,
    'BA' :  56,
    'La' :  57,
    'LA' :  57,
    'Ce' :  58,
    'CE' :  58,
    'Pr' :  59,
    'PR' :  59,
    'Nd' :  60,
    'ND' :  60,
    'Pm' :  61,
    'PM' :  61,
    'Sm' :  62,
    'SM' :  62,
    'Eu' :  63,
    'EU' :  63,
    'Gd' :  64,
    'GD' :  64,
    'Tb' :  65,
    'TB' :  65,
    'Dy' :  66,
    'DY' :  66,
    'Ho' :  67,
    'HO' :  67,
    'Er' :  68,
    'ER' :  68,
    'Tm' :  69,
    'TM' :  69,
    'Yb' :  70,
    'YB' :  70,
    'Lu' :  71,
    'LU' :  71,
    'Hf' :  72,
    'HF' :  72,
    'Ta' :  73,
    'TA' :  73,
    'W'  :  74,
    'Re' :  75,
    'RE' :  75,
    'Os' :  76,
    'OS' :  76,
    'Ir' :  77,
    'IR' :  77,
    'Pt' :  78,
    'PT' :  78,
    'Au' :  79,
    'AU' :  79,
    'Hg' :  80,
    'HG' :  80,
    'Tl' :  81,
    'TL' :  81,
    'Pb' :  82,
    'PB' :  82,
    'Bi' :  83,
    'BI' :  83,
    'Po' :  84,
    'PO' :  84,
    'At' :  85,
    'AT' :  85,
    'Rn' :  86,
    'RN' :  86,
    'Fr' :  87,
    'FR' :  87,
    'Ra' :  88,
    'RA' :  88,
    'Ac' :  89,
    'AC' :  89,
    'Th' :  90,
    'TH' :  90,
    'Pa' :  91,
    'PA' :  91,
    'U'  :  92,
    'Np' :  93,
    'NP' :  93,
    'Pu' :  94,
    'PU' :  94,
    'Am' :  95,
    'AM' :  95,
    'Cm' :  96,
    'CM' :  96,
    'Bk' :  97,
    'BK' :  97,
    'Cf' :  98,
    'CF' :  98,
    'Es' :  99,
    'ES' :  99,
    'Fm' : 100,
    'FM' : 100,
    'Md' : 101,
    'MD' : 101,
    'No' : 102,
    'NO' : 102,
    'Lr' : 103,
    'LR' : 103,
    'Rf' : 104,
    'RF' : 104,
    'Db' : 105,
    'DB' : 105,
    'Sg' : 106,
    'SG' : 106,
    'Bh' : 107,
    'BH' : 107,
    'Hs' : 108,
    'HS' : 108,
    'Mt' : 109,
    'MT' : 109
    }

implicit_valence = {
    'H'  :  {0:1,1:0},
    'C'  :  {0:4,1:3,2:2,3:1,4:0},
    'N'  :  {0:3,1:2,2:1,3:0},
    'O'  :  {0:2,1:1,2:0},
    'F'  :  {0:1,1:0},
    'Cl' :  {0:1,1:0},
    'CL' :  {0:1,1:0},
    'Br' :  {0:1,1:0},
    'BR' :  {0:1,1:0},
    'I'  :  {0:1,1:0},
    'S'  :  {0:2,1:2,2:0,3:1,4:0,5:1,6:0}, # ambiguity?
    'K'  :  {0:0,1:0},     # as drawn
    'Cu' :  {0:0,1:0,2:0}, # as drawn
    'CU' :  {0:0,1:0,2:0}, # as drawn
    'Zn' :  {0:0,1:0,2:0}, # as drawn
    'ZN' :  {0:0,1:0,2:0}, # as drawn
    'Mg' :  {0:1,1:0},
    'MG' :  {0:1,1:0},
    'Ca' :  {0:1,1:0},
    'CA' :  {0:1,1:0},
    'P'  :  {6:0},
    }
