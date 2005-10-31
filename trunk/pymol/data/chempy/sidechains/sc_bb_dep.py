# pymol -c generate2.py

from chempy import io
from glob import glob
from copy import deepcopy

# backbone-independent rotamers

lines = io.lst.fromFile("bbdep02.May.lib")

chi = { 'CYS' :
        { '1': ('N'  , 'CA' , 'CB' , 'SG' )   },
        'ASP' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'OD1'), },
        'GLU' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD' ), 
          '3': ('CB' , 'CG' , 'CD' , 'OE1'), },
        'PHE' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD1'), },
        'HIS' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'ND1'), },
        'ILE' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG1'),
          '2': ('CA' , 'CB' , 'CG1', 'CD1+CD'), }, 
        'LYS' :
        { '1': ('N'  , 'CA' , 'CB'  ,'CG' ),
          '2': ('CA' , 'CB' , 'CG'  ,'CD' ),
          '3': ('CB' , 'CG' , 'CD'  ,'CE' ), 
          '4': ('CG' , 'CD' , 'CE'  ,'NZ' ), },
        'LEU' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD1'), }, 
        'MET' :
        { '1': ('N'  , 'CA' , 'CB'  ,'CG' ),
          '2': ('CA' , 'CB' , 'CG'  ,'SD' ),
          '3': ('CB' , 'CG' , 'SD'  ,'CE' ), },
        'ASN' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'OD1'), },
        'PRO' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD' ), }, 
        'GLN' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD' ), 
          '3': ('CB' , 'CG' , 'CD' , 'OE1'), },
        'ARG' :
        { '1': ('N'  , 'CA' , 'CB'  ,'CG' ),
          '2': ('CA' , 'CB' , 'CG'  ,'CD' ),
          '3': ('CB' , 'CG' , 'CD'  ,'NE' ), 
          '4': ('CG' , 'CD' , 'NE'  ,'CZ' ), },
        'SER' :
        { '1': ('N'  , 'CA' , 'CB' , 'OG' ), },
        'THR' :
        { '1': ('N'  , 'CA' , 'CB' , 'OG1'), },
        'VAL' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG1'), },
        'TRP' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD1'), },
        'TYR' :
        { '1': ('N'  , 'CA' , 'CB' , 'CG' ),
          '2': ('CA' , 'CB' , 'CG' , 'CD1'), },
        }


total = {}
output = {}

# first, total the number of rotamers for each residue

for line in lines:
    field = line.split()

    if len(field)>6:
        resn = field[0]
        phi = int(field[1])
        psi = int(field[2])
        key = (resn, phi, psi)
        total[key] = total.get(key,0) + float(field[8])

# now build library, including frequency value

dd  = {}

for line in lines:
    field = line.split()
    len_field = len(field)
    if len_field>11:
            resn = field[0]
            phi = int(field[1])
            psi = int(field[2])
            if (((psi==(60*(psi/60))) and (phi==(60*(phi/60)))) or 
                ((psi==(20*(psi/20))) and (phi==(20*(phi/20))) and (int(field[3])>0)) or
                ((psi==(10*(psi/10))) and (phi==(10*(phi/10))) and (int(field[3])>=250))):
                if not dd.has_key((phi,psi)):
#                    print phi,psi
                    dd[(phi,psi)]=1
                key = (resn, phi, psi)
                list = output.get(key,[])
                fsum = float(field[8])
                dict = {}
                if total[key]>0.0:
                    freq = fsum/total[key]
                else:
                    freq = 0.0
                if freq>0.01:
                    dict['FREQ'] = freq
                    len_chi = len(chi[resn].keys())
                    if len_chi>0:
                        chi1 = float(field[9])
                        dict[chi[resn]['1']] = chi1
                    if len_chi>1:
                        chi2 = float(field[10])
                        dict[chi[resn]['2']] = chi2            
                    if len_chi>2:
                        chi3 = float(field[11])
                        dict[chi[resn]['3']] = chi3                        
                    if len_chi>3:
                        chi4 = float(field[12])
                        dict[chi[resn]['4']] = chi4
                    list.append((freq,dict))
                    output[key] = list

# sort by priority

for key in output.keys():
    list = output[key]
    list.sort()
    list.reverse()
    output[key] = map(lambda x:x[1],list)

#for key in output.keys():
#    
#    if key[0] == 'HIS':
#        output[( 'HIE', ) + key[1:]] = output[key]
#        output[( 'HID', ) + key[1:]] = output[key]
#        output[( 'HIP', ) + key[1:]] = output[key]        
#    elif key[0] == 'CYS':
#        output[( 'CYX', ) + key[1:]] = output[key]                

io.pkl.toFile(output,"sc_bb_dep.pkl")


