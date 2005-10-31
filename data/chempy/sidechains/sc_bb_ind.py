# pymol -c generate2.py

from chempy import io
from glob import glob
from copy import deepcopy

# backbone-independent rotamers

lines = io.lst.fromFile("bbind02.May.lib")

# skip to the data section

while lines[0][0:3]!='Res':
    lines.pop(0)
lines.pop(0)

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
        total[resn] = total.get(resn,0) + int(field[6])

# now build library, including frequency value

for line in lines:
    field = line.split()
    len_field = len(field)
    if len_field>11:
        resn = field[0]
        list = output.get(resn,[])
        count = float(field[6])
        dict = {}
        freq = count/total[resn]
        dict['FREQ'] = freq
        if len_field>11:
            chi1 = float(field[11])
            dict[chi[resn]['1']] = chi1
        if len_field>13:
            chi2 = float(field[13])
            dict[chi[resn]['2']] = chi2            
        if len_field>15:
            chi3 = float(field[15])
            dict[chi[resn]['3']] = chi3                        
        if len_field>17:
            chi4 = float(field[17])
            dict[chi[resn]['4']] = chi4
        list.append((freq,dict))
        output[resn] = list

# sort by priority

for resn in output:
    list = output[resn]
    list.sort()
    list.reverse()
    output[resn] = map(lambda x:x[1],list)

# common aliases

#output['HIE'] = output['HIS']
#output['HID'] = output['HIS']
#output['HIP'] = output['HIS']
#output['CYX'] = output['CYS']
    
io.pkl.toFile(output,"sc_bb_ind.pkl")


