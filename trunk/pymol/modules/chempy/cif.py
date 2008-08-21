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

from chempy import io, Atom, Bond
from chempy.models import Indexed

# OMG what a horrid file format.

single_quote_re = re.compile(r"'[^']*[']*'") # doesn't yet handle ESC
double_quote_re = re.compile(r'"[^"]*["]*"') # ditto
bracket_quote_re = re.compile(r'\[[^\]]*\]') # ditto

clean_float_re = re.compile(r'[^0-9+\-.eE].*')
clean_int_re = re.compile(r'[^0-9+\-].*')

class CIFRec:

    def get_quoted_value(self):
        if self.line[0:1] == "'":
            mo = single_quote_re.match(self.line)
        elif self.line[0:1] == '"':
            mo = double_quote_re.match(self.line)
        elif self.line[0:1] == '[':
            mo = bracket_quote_re.match(self.line)
        else:
            mo = None
        if mo != None:
            result = self.line[:mo.end()]
            self.line = self.line[mo.end():]
            if len(string.strip(self.line[0:1])): # followed by non-whitespace...
                self.line = result[-1:] + self.line
                result = result[:-1] + self.get_quoted_value()
            self.check_line()
        else:
            result = None # shouldn't happen
        return result
    
    def get_delimited_value(self):
        result = ''
        while 1:
            if self.line[0:1] == ';':
                self.line = self.line[1:]
                self.check_line()
                break
            result = result + self.line
            self.next_line()
        return result

    def get_next_word(self):
        result = None
        mo = re.search("\s+", self.line)
        if mo == None:
            result = self.line
            self.next_line()
        else:
            result = self.line[:mo.start()]
            self.line = self.line[mo.end():]
            self.check_line()
        return result
    
    def trim_leading_whitespace(self):
        while self.line != None:
            mo = re.match("\s",self.line)
            if mo != None:
                self.line = self.line[mo.end():]
                self.check_line()
            else:
                break

    def get_next_value(self):
        self.trim_leading_whitespace()
        if self.line != None:
            if self.line[0:1] == '_':
                return None
            elif self.line[0:1] == ';':
                self.line = self.line[1:]
                self.check_line()
                return self.get_delimited_value()
            elif string.lower(self.line[0:5]) in ['loop_','data_','save_']:
                return None
            elif string.lower(self.line[0:6]) == 'GLOBAL':
                return None
            elif self.line[0:1] in [ "'", '"']: #  '[' ]: bracket quote fubar with PDB's mmCIF data???
                return self.get_quoted_value()
            else:
                return self.get_next_word()
        return None
        
    def check_line(self):
        if self.line != None:
            while not len(string.strip(self.line)):
                self.next_line()
                if self.line == None:
                    break

    def next_line(self):
        self.line = None
        while len(self.list):
            self.line = self.list.pop()
            # nuke hash comments 
            hash = string.find(self.line,'#')
            if hash>=0:
                self.line = self.line[0:hash]
            # and only return non-blank lines
            if len(string.strip(self.line)):
                break
#        print "next_line:", self.line,
        
    def parse_loop_body(self,fields):
        len_fields = len(fields)
        records = []
        record = []
        cnt = len_fields
        while self.line != None: 
            value = self.get_next_value()
            if value == None:
                break
            else:
                record.append(value)
                cnt = cnt - 1
                if not cnt:
                    cnt = len_fields
                    if len(record) == len_fields:
                        records.append(record)
                    record = []
#                print "loop_read [%s]=[%s]"%(fields[0],value)
        if len(record) == len_fields:
            records.append(record)
        self.loops.append( (fields,records) )
        
    def parse_loop(self):
#        print "parsing loop..."
        fields = []
        while self.line != None: 
            if string.lower(self.line[0:5])=='loop_':
                break
            elif self.line[0:1]=='_':
                fields.append(string.lower(string.strip(self.line)))
            else:
                self.parse_loop_body(fields)
                break
            self.next_line()
                
    def parse_name_value(self):
        data_name = self.get_next_word()
        data_value = self.get_next_value()
        self.key_value[data_name] = data_value
#        print "data_read [%s]=[%s]"%(data_name,data_value)
        
    def parse_normal(self):
        while self.line != None:
            if self.line[0:1] == '_': # data name
                self.parse_name_value()
            elif string.lower(self.line[0:5])=='loop_':
                print "entering loop",self.line
                self.line = self.line[5:]
                self.check_line()
                self.parse_loop()
                print "exiting loop",self.line                
            else: # shouldn't happen
                print "unhandled: [%s]"%self.line
                self.next_line()

    def to_float(self, key):
        value = self.key_value[key]
        value = clean_float_re.sub("",value)
        return float(value)

    def to_str(self, key):
        value = self.key_value[key]
        if value[0:1]=="'" and value[-1:] =="'":
            value = value[1:-1]
        if value[0:1]=='"' and value[-1:] =='"':
            value = value[1:-1]
        return value
    
    def read_symmetry(self):
        kv = self.key_value

        if ( kv.has_key("_cell_length_a") and
             kv.has_key("_cell_length_b") and
             kv.has_key("_cell_length_c") and
             kv.has_key("_cell_angle_alpha") and
             kv.has_key("_cell_angle_beta") and
             kv.has_key("_cell_angle_gamma")):
            self.model.cell = [
                self.to_float("_cell_length_a"),
                self.to_float("_cell_length_b"),
                self.to_float("_cell_length_c"),
                self.to_float("_cell_angle_alpha"),
                self.to_float("_cell_angle_beta"),
                self.to_float("_cell_angle_gamma")]
#        if hasattr(self.model,'cell'):
#            print self.model.cell
        if kv.has_key("_symmetry_space_group_name_h-m"):
            self.model.spacegroup = self.to_str('_symmetry_space_group_name_h-m')
#        if hasattr(self.model,'spacegroup'):
#            print self.model.spacegroup

    def index_to_int(self, index, value):
        result = clean_int_re.sub("",value[index])
        return int(result)

    def index_to_float(self, index, value):
        result = clean_float_re.sub("",value[index])
        if len(result):
            return float(result)
        else:
            return 0.0
        

    def index_to_str(self, index, value):
        return value[index]
                
    def read_chem_comp_atom_model_cartn(self,fields,field_dict,values):
        cartn_x = field_dict['_chem_comp_atom_model_cartn_x']
        cartn_y = field_dict['_chem_comp_atom_model_cartn_y']
        cartn_z = field_dict['_chem_comp_atom_model_cartn_z']
        print cartn_x,cartn_y,cartn_z
        name = field_dict.get('_chem_comp_atom_atom_id',None)
        symbol = field_dict.get('_chem_comp_atom_type_symbol',None)
        resn = field_dict.get('_chem_comp_atom.comp_id',None)
        partial_charge = field_dict.get('_chem_comp_atom_partial_charge',None)
        formal_charge = field_dict.get('chem_comp_atom_charge',None)
        str_fields = []
        if symbol != None: str_fields.append( ('symbol',symbol) )
        if name != None: str_fields.append( ('name',name) )
        if resn != None: str_fields.append( ('resn',resn) )
        float_fields = []
        if partial_charge != None: float_fields.append( ('partial_charge',partial_charge) )
        int_fields = []
        if formal_charge != None: int_fields.append( ('formal_charge',formal_charge) ) 
        for value in values:
            atom = Atom()
            atom.coord = [
                self.index_to_float(cartn_x,value),
                self.index_to_float(cartn_y,value),
                self.index_to_float(cartn_z,value)]
            self.model.atom.append(atom)
            for field in str_fields:
                setattr(atom,field[0],value[field[1]])
            for field in float_fields:
                setattr(atom,field[0],self.index_to_float(field[1],value))
            for field in int_fields:
                setattr(atom,field[0],self.index_to_int(field[1],value))

    def read_chem_comp_atom(self):
        self.atom_site_label_index = {}
        for loop in self.loops:
            (fields, field_dict, values) = loop
            if (field_dict.has_key("_chem_comp_atom_model_cartn_x") and
                field_dict.has_key("_chem_comp_atom_model_cartn_y") and
                field_dict.has_key("_chem_comp_atom_model_cartn_z")): 
                self.read_chem_comp_atom_model_cartn(fields,field_dict,values)

    def read_atom_site_fract(self,fields,field_dict,values):
        self.model.fractional = 1

        fract_x = field_dict['_atom_site_fract_x']
        fract_y = field_dict['_atom_site_fract_y']
        fract_z = field_dict['_atom_site_fract_z']
        symbol = field_dict.get('_atom_site_type_symbol',None)
        name = field_dict.get('_atom_site_label',None)
        u = field_dict.get('_atom_site_u_iso_or_equiv',None)                
        str_fields = []
        if symbol != None: str_fields.append( ('symbol',symbol) )
        if name != None: str_fields.append( ('name',name) )
        float_fields = []
        if u != None: float_fields.append( ('u', u))
        int_fields = []
        for value in values:
            atom = Atom()
            atom.coord = [
                self.index_to_float(fract_x,value),
                self.index_to_float(fract_y,value),
                self.index_to_float(fract_z,value)]
            self.model.atom.append(atom)
            for field in str_fields:
                setattr(atom,field[0],value[field[1]])
            for field in float_fields:
                setattr(atom,field[0],self.index_to_float(field[1],value))
            for field in int_fields:
                setattr(atom,field[0],self.index_to_int(field[1],value))

    def read_atom_site_cartn(self,fields,field_dict,values):
        cartn_x = field_dict['_atom_site_cartn_x']
        cartn_y = field_dict['_atom_site_cartn_y']
        cartn_z = field_dict['_atom_site_cartn_z']
        group_pdb = field_dict.get('_atom_site_group_pdb',None)
        symbol = field_dict.get('_atom_site_type_symbol',None)
        name = field_dict.get('_atom_site_label_atom_id',None)
        resn = field_dict.get('_atom_site_label_comp_id',None)
        resi = field_dict.get('_atom_site_label_seq_id',None)
        chain = field_dict.get('_atom_site_label_asym_id',None)
        ins_code = field_dict.get('_atom_site_pdbx_pdb_ins_code',None)
        # use auth fields preferentially, if provided
        auth_resn = field_dict.get('_atom_site_auth_comp_id',None)
        auth_resi = field_dict.get('_atom_site_auth_seq_id',None)
        auth_name = field_dict.get('_atom_site_auth_atom_id',None)
        auth_chain = field_dict.get('_atom_site_auth_asym_id',None)        
        if auth_resn != None: resn = auth_resn
        if auth_resi != None: resi = auth_resi
        if auth_name != None: name = auth_name
        if auth_chain != None: chain = auth_chain
        b = field_dict.get('_atom_site_b_iso_or_equiv',None)        
        q = field_dict.get('_atom_site_occupancy',None)
        ID = field_dict.get('_atom_site_id',None)
        str_fields = []
        if symbol != None: str_fields.append( ('symbol',symbol) )
        if name != None: str_fields.append( ('name',name) )
        if resn != None: str_fields.append( ('resn',resn) )
        if resi != None: str_fields.append( ('resi',resi) )
        if chain != None: str_fields.append( ('chain',chain) )        
        if ins_code != None: str_fields.append( ('ins_code',ins_code) )
        float_fields = []
        if q != None: float_fields.append( ('q',q) )
        if b != None: float_fields.append( ('b',b) )
        int_fields = []
        if ID != None: int_fields.append( ('id',ID) )                                

        for value in values:
            atom = Atom()
            atom.coord = [
                self.index_to_float(cartn_x,value),
                self.index_to_float(cartn_y,value),
                self.index_to_float(cartn_z,value)]
            self.model.atom.append(atom)
            if group_pdb != None:
                if value[group_pdb] == 'ATOM':
                    atom.hetatm = 0
                else:
                    atom.hetatm = 1
            for field in str_fields:
                setattr(atom,field[0],value[field[1]])
            for field in float_fields:
                setattr(atom,field[0],self.index_to_float(field[1],value))
            for field in int_fields:
                setattr(atom,field[0],self.index_to_int(field[1],value))
                
#        for a in self.model.atom:
#            print a.coord
        
    def read_atom_site_aniso(self,fields,field_dict,values):
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        label = field_dict['_atom_site_aniso_label']
        u11 = field_dict['_atom_site_aniso_u_11']
        u22 = field_dict['_atom_site_aniso_u_22']
        u33 = field_dict['_atom_site_aniso_u_33']
        u12 = field_dict['_atom_site_aniso_u_12']
        u13 = field_dict['_atom_site_aniso_u_13']
        u23 = field_dict['_atom_site_aniso_u_23']
        for value in values:
            atom = name_dict[self.index_to_str(label,value)]
            self.model.atom[atom].u_aniso = [
                self.index_to_float(u11,value),
                self.index_to_float(u22,value),
                self.index_to_float(u33,value),
                self.index_to_float(u12,value),
                self.index_to_float(u13,value),
                self.index_to_float(u23,value)
                ]

    def read_atom_site(self):
        self.atom_site_label_index = {}
        for loop in self.loops:
            (fields, field_dict, values) = loop
            if (field_dict.has_key("_atom_site_fract_x") and
                field_dict.has_key("_atom_site_fract_y") and
                field_dict.has_key("_atom_site_fract_z")): # fractional coords
                self.read_atom_site_fract(fields,field_dict,values)
            elif (field_dict.has_key("_atom_site_cartn_x") and
                  field_dict.has_key("_atom_site_cartn_y") and
                  field_dict.has_key("_atom_site_cartn_z")): # cartesian coords
                self.read_atom_site_cartn(fields,field_dict,values)
            elif (field_dict.has_key("_atom_site_aniso_label") and
                  field_dict.has_key("_atom_site_aniso_u_11") and
                  field_dict.has_key("_atom_site_aniso_u_22") and
                  field_dict.has_key("_atom_site_aniso_u_33") and
                  field_dict.has_key("_atom_site_aniso_u_12") and
                  field_dict.has_key("_atom_site_aniso_u_13") and
                  field_dict.has_key("_atom_site_aniso_u_23")): # anisotropics
                self.read_atom_site_aniso(fields,field_dict,values)

    def read_geom_bond_atom_site_labels(self,fields,field_dict,values):
        # create index of atom name
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        label_1 = field_dict['_geom_bond_atom_site_label_1']
        label_2 = field_dict['_geom_bond_atom_site_label_2']
        for value in values:
            bond = Bond()
            bond.index = [
                name_dict[self.index_to_str(label_1,value)],
                name_dict[self.index_to_str(label_2,value)]]
            bond.order = 1
            self.model.bond.append(bond)

    def read_geom_bond(self):
        self.atom_site_label_index = {}
        for loop in self.loops:
            (fields, field_dict, values) = loop
            if (field_dict.has_key("_geom_bond_atom_site_label_1") and
                field_dict.has_key("_geom_bond_atom_site_label_2")):
                self.read_geom_bond_atom_site_labels(fields,field_dict,values)
    
    def read_chem_comp_bond_atom_ids(self,fields,field_dict,values):
        order_table = { 'sing' : 1, 'doub' : 2, 'trip' :3, 'delo': 4 }
        # create index of atom name
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        label_1 = field_dict['_chem_comp_bond_atom_id_1']
        label_2 = field_dict['_chem_comp_bond_atom_id_2']
        order = field_dict.get('_chem_comp_bond_value_order',None)
        for value in values:
            bond = Bond()
            bond.index = [
                name_dict[self.index_to_str(label_1,value)],
                name_dict[self.index_to_str(label_2,value)]]
            if order != None:
                order_string = string.lower(self.index_to_str(order,value))
                bond.order = order_table.get(order_string[0:4],1)
            else:
                bond.order = 1
            self.model.bond.append(bond)

    def read_chem_comp_bond(self):
        self.atom_site_label_index = {}
        for loop in self.loops:
            (fields, field_dict, values) = loop
            if (field_dict.has_key("_chem_comp_bond_atom_id_1") and
                field_dict.has_key("_chem_comp_bond_atom_id_2")):
                self.read_chem_comp_bond_atom_ids(fields,field_dict,values)

    def create_field_dicts(self):
        new_loops = []
        for loop in self.loops:
            (fields, values) = loop
            field_dict = {}
            cnt = 0
            for field in fields:
                field_dict[field] = cnt
                cnt = cnt + 1
            new_loops.append( (fields,field_dict,values) )
        self.loops = new_loops
                
    def convert_dot_to_underscore(self):
        for key in self.key_value.keys():
            if string.find(key,".")>=0:
                new_key = string.replace(key,".","_")
            else:
                new_key = key
            new_key = string.lower(new_key)
            if new_key != key:
                self.key_value[new_key] = self.key_value[key]
                del self.key_value[key]
        new_loops = []
        for loop in self.loops:
            (fields, values) = loop
            fields = map(lambda x:string.replace(x,".","_"),fields)
            fields = map(lambda x:string.lower(x),fields)            
            new_loops.append( (fields,values) )
        self.loops = new_loops
        
    def __init__(self,cif_list):
        cif_list.reverse()
        self.list = cif_list
        data_line = self.list.pop()
        self.loops = []
        self.key_value = {}
        self.data_name = string.strip(data_line[5:])
        self.next_line()
        self.parse_normal()
        print ' CIF: For data block "%s"...'%self.data_name
        print " CIF: Read %d key/value pair(s)."%len(self.key_value)
        print " CIF: Read %d table(s)."%len(self.loops)

        self.convert_dot_to_underscore()
        self.create_field_dicts()

        # now build the molecule record
        self.model = Indexed()

        # by default, indicate that we want PyMOL to automatically
        # detect bonds based on coordinates

        self.model.connect_mode = 3
        
        self.read_symmetry()
        self.read_atom_site()
        self.read_chem_comp_atom()

        self.read_geom_bond()
        self.read_chem_comp_bond()
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
        data_seen = 0
        while 1:
            if self.input_line == None:
                s = self.file.readline()
            else:
                s = self.input_line
                self.input_line = None
            if not s: # end of file
                if len(cur)>0:
                    return CIFRec(cur)
                else:
                    return None
            elif string.lower(s[0:5]) == r'data_': # signals start of new record
                data_seen = 1
                if len(cur)>0:
                    self.input_line = s # save for next time
                    return CIFRec(cur)
                elif data_seen:
                    cur.append(s)
            elif data_seen:
                cur.append(s)
    def close(self):
        self.file.close()
        
    
