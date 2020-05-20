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

import re

from chempy import Atom, Bond
from chempy.models import Indexed

# OMG what a horrid file format.

# set asym_id as chain AND segi
ASYM_ID_AS_SEGI = True

# matches any valid CIF token
token_re = re.compile(r'(?:\s*('
    r"'.*?'(?=\s)|"
    r'".*?"(?=\s)|'
    r'^;.*?^;|'
    r'#.*?$|'
    r'\S+'
    r'))', re.I | re.DOTALL | re.MULTILINE)

# matches numbers in parentheses
floatuncert_sub = re.compile(r'\([0-9]+\)').sub

def unquote(s):
    '''
    Remove quotes from a CIF token
    '''
    s0 = s[0]
    if s0 in ('"', "'"):
        return s[1:-1]
    if s0 == ';':
        return s[1:-2]
    if s0 == '[':
        raise ValueError('found reserved "[" character')
    return s

def scifloat(s):
    '''
    Convert string to float. Accepts numbers which match the CIF regex
    for floats, which may contain an uncertenty in parentheses
    '''
    # strip uncertenty notation from string: 123(4)E02 -> 123E02
    return float(floatuncert_sub('', s))

class ciftokeniter(object):
    '''
    Tokenize a CIF string. Skip comments, preserve quotes. Provides
    sub-iterators for loop keys and data.
    '''
    def __init__(self, starstr):
        self._iter = token_re.finditer(starstr)
        self._prev = []
    def __iter__(self):
        return self
    def next(self):
        if self._prev:
            return self._prev.pop()
        s = next(self._iter).group(1)
        if s[0] == '#':
            return next(self)
        return s
    __next__ = next
    def loopdataiter(self):
        for s in self:
            if s[0] == '_' or s[:5].lower() in ('loop_', 'data_', 'save_'):
                self._prev.append(s)
                break
            yield s
    def loopkeysiter(self):
        for s in self:
            if s[0] != '_':
                self._prev.append(s)
                break
            yield s

def parse_cif(cifstr):
    '''
    Parse a CIF string and return an iterator over CIFData records.
    '''
    current_data = current_block = None

    token_it = ciftokeniter(cifstr)

    for s in token_it:
        s_lower = s.lower()
        if s[0] == '_':
            key = s_lower.replace('.', '_')
            current_block.key_value[key] = next(token_it)
        elif s_lower == 'loop_':
            loop = CIFLoop()
            current_block.loops.append(loop)
            for i, key in enumerate(token_it.loopkeysiter()):
                key = key.lower().replace('.', '_')
                loop.keys[key] = i
            ncols = len(loop.keys)
            for i, value in enumerate(token_it.loopdataiter()):
                if i % ncols == 0:
                    row = []
                    loop.rows.append(row)
                row.append(value)
        elif s_lower[:5] == 'data_':
            if current_data is not None:
                yield current_data
            current_block = current_data = CIFData(s[5:])
        elif s_lower[:5] == 'save_':
            if len(s) > 5:
                current_block = CIFData(s[5:])
                current_data.saveframes.append(current_block)
            else:
                current_block = current_data
        else:
            raise ValueError(s)

    if current_data is not None:
        yield current_data

def _row_get(row, i, d, cast):
    if i < 0:
        return d
    v = row[i]
    if v in ('.', '?'):
        return d
    return cast(v)

class CIFLoop:
    '''
    CIF loop (table)
    '''
    def __init__(self):
        self.keys = {}
        self.rows = []

    def get_col_idx_opt(self, *names):
        '''
        Get column index for first found name in names, or -1
        '''
        for name in names:
            idx = self.keys.get(name.replace('.', '_'))
            if idx is not None:
                return idx
        return -1

    def get_col_idx(self, *names):
        '''
        Get column index for first found name in names, or raise KeyError
        '''
        idx = self.get_col_idx_opt(*names)
        if idx == -1:
            raise KeyError
        return idx

class CIFData:
    '''
    CIF data
    '''
    def __init__(self, name):
        self.name = name
        self.loops = []
        self.key_value = {}
        self.saveframes = {}

    def __repr__(self):
        return '<%s:%s #kv=%d #loops=%d>' % (type(self).__name__,
                self.name, len(self.key_value), len(self.loops))

    def _get(self, key, d, cast):
        v = self.key_value[key]
        if v in ('.', '?'):
            return d
        return cast(v)

    def to_float(self, key):
        return self._get(key, 0.0, scifloat)

    def to_str(self, key):
        return self._get(key, '', unquote)

    def index_to_int(self, index, value):
        return _row_get(value, index, 0, int)

    def index_to_float(self, index, value):
        return _row_get(value, index, 0.0, scifloat)

    def index_to_str(self, index, value):
        return _row_get(value, index, '', unquote)

class CIFRec(CIFData):
    '''
    CIF record
    '''
    def __init__(self, datablock):
        self.loops = datablock.loops
        self.key_value = datablock.key_value
        self.name = datablock.name

        # now build the molecule record
        self.model = Indexed()
        self.model.molecule.title = datablock.name

        # coordinates for state 2-N
        self.extra_coords = []

        # by default, indicate that we want PyMOL to automatically
        # detect bonds based on coordinates
        self.model.connect_mode = 3

        self.read_symmetry()

        if self.read_atom_site():
            self.read_geom_bond()
            self.read_struct_conn()
            self.read_ss()
        elif self.read_chem_comp_atom():
            self.read_chem_comp_bond()

    def read_symmetry(self):
        try:
            self.model.cell = [
                self.to_float("_cell_length_a"),
                self.to_float("_cell_length_b"),
                self.to_float("_cell_length_c"),
                self.to_float("_cell_angle_alpha"),
                self.to_float("_cell_angle_beta"),
                self.to_float("_cell_angle_gamma")]
        except (KeyError, ValueError):
            return False

        try:
            self.model.spacegroup = self.to_str('_symmetry_space_group_name_h-m')
        except KeyError:
            pass

        return True

    def read_chem_comp_atom_model_cartn(self,fields,field_dict,values):
        try:
            cartn_x = field_dict['_chem_comp_atom_model_cartn_x']
            cartn_y = field_dict['_chem_comp_atom_model_cartn_y']
            cartn_z = field_dict['_chem_comp_atom_model_cartn_z']
        except KeyError:
            return False
        name = field_dict.get('_chem_comp_atom_atom_id',None)
        symbol = field_dict.get('_chem_comp_atom_type_symbol',None)
        resn = field_dict.get('_chem_comp_atom_comp_id',None)
        partial_charge = field_dict.get('_chem_comp_atom_partial_charge',None)
        formal_charge = field_dict.get('chem_comp_atom_charge',None)
        str_fields = []
        if symbol is not None: str_fields.append( ('symbol',symbol) )
        if name is not None: str_fields.append( ('name',name) )
        if resn is not None: str_fields.append( ('resn',resn) )
        float_fields = []
        if partial_charge is not None: float_fields.append( ('partial_charge',partial_charge) )
        int_fields = []
        if formal_charge is not None: int_fields.append( ('formal_charge',formal_charge) )
        for value in values:
            atom = Atom()
            atom.coord = [
                self.index_to_float(cartn_x,value),
                self.index_to_float(cartn_y,value),
                self.index_to_float(cartn_z,value)]
            self.model.atom.append(atom)
            for field in str_fields:
                setattr(atom,field[0],self.index_to_str(field[1],value))
            for field in float_fields:
                setattr(atom,field[0],self.index_to_float(field[1],value))
            for field in int_fields:
                setattr(atom,field[0],self.index_to_int(field[1],value))
        return True

    def read_chem_comp_atom(self):
        for loop in self.loops:
            if self.read_chem_comp_atom_model_cartn(None, loop.keys, loop.rows):
                return True
        return False

    def read_atom_site_fract(self,fields,field_dict,values):
        try:
            fract_x = field_dict['_atom_site_fract_x']
            fract_y = field_dict['_atom_site_fract_y']
            fract_z = field_dict['_atom_site_fract_z']
        except KeyError:
            return False
        self.model.fractional = 1
        symbol = field_dict.get('_atom_site_type_symbol',None)
        name = field_dict.get('_atom_site_label',None)
        if name is None:
            name = field_dict.get('_atom_site_id',None)
        u = field_dict.get('_atom_site_u_iso_or_equiv',None)
        str_fields = []
        if symbol is not None: str_fields.append( ('symbol',symbol) )
        if name is not None: str_fields.append( ('name',name) )
        float_fields = []
        if u is not None: float_fields.append( ('u', u))
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
        return True

    def read_atom_site_cartn(self,fields,field_dict,values):
        try:
            cartn_x = field_dict['_atom_site_cartn_x']
            cartn_y = field_dict['_atom_site_cartn_y']
            cartn_z = field_dict['_atom_site_cartn_z']
        except KeyError as e:
            return False
        group_pdb = field_dict.get('_atom_site_group_pdb',None)
        symbol = field_dict.get('_atom_site_type_symbol',None)
        name = field_dict.get('_atom_site_label_atom_id',None)
        resn = field_dict.get('_atom_site_label_comp_id',None)
        resi = field_dict.get('_atom_site_label_seq_id',None)
        chain = field_dict.get('_atom_site_label_asym_id',None)
        ins_code = field_dict.get('_atom_site_pdbx_pdb_ins_code',None)
        alt = field_dict.get('_atom_site_label_alt_id', None)
        model_num = field_dict.get('_atom_site_pdbx_pdb_model_num', None)
        # use auth fields preferentially, if provided
        auth_resn = field_dict.get('_atom_site_auth_comp_id',None)
        auth_resi = field_dict.get('_atom_site_auth_seq_id',None)
        auth_name = field_dict.get('_atom_site_auth_atom_id',None)
        auth_chain = field_dict.get('_atom_site_auth_asym_id',None)
        if auth_resn is not None: resn = auth_resn
        if auth_resi is not None: resi = auth_resi
        if auth_name is not None: name = auth_name
        if auth_chain is not None: chain = auth_chain
        b = field_dict.get('_atom_site_b_iso_or_equiv',None)
        q = field_dict.get('_atom_site_occupancy',None)
        ID = field_dict.get('_atom_site_id',None)
        str_fields = []
        if symbol is not None: str_fields.append( ('symbol',symbol) )
        if name is not None: str_fields.append( ('name',name) )
        if resn is not None: str_fields.append( ('resn',resn) )
        if resi is not None: str_fields.append( ('resi',resi) )
        if chain is not None:
            str_fields.append( ('chain',chain) )
            if ASYM_ID_AS_SEGI:
                str_fields.append( ('segi',chain) )
        if alt is not None: str_fields.append( ('alt',alt) )
        if ins_code is not None: str_fields.append( ('ins_code',ins_code) )
        float_fields = []
        if q is not None: float_fields.append( ('q',q) )
        if b is not None: float_fields.append( ('b',b) )
        int_fields = []
        if ID is not None: int_fields.append( ('id',ID) )

        first_model_num = self.index_to_int(model_num, values[0])

        for value in values:
            coord = [
                self.index_to_float(cartn_x,value),
                self.index_to_float(cartn_y,value),
                self.index_to_float(cartn_z,value)]

            if model_num is not None:
                v = self.index_to_int(model_num, value)
                if v != first_model_num:
                    self.extra_coords.extend(coord)
                    continue

            atom = Atom()
            atom.coord = coord
            self.model.atom.append(atom)
            if group_pdb is not None:
                if value[group_pdb] == 'ATOM':
                    atom.hetatm = 0
                else:
                    atom.hetatm = 1
            for field in str_fields:
                setattr(atom,field[0],self.index_to_str(field[1],value))
            for field in float_fields:
                setattr(atom,field[0],self.index_to_float(field[1],value))
            for field in int_fields:
                setattr(atom,field[0],self.index_to_int(field[1],value))
        return True

    def read_atom_site_aniso(self,fields,field_dict,values):
        try:
            label = field_dict['_atom_site_aniso_label']
            u11 = field_dict['_atom_site_aniso_u_11']
            u22 = field_dict['_atom_site_aniso_u_22']
            u33 = field_dict['_atom_site_aniso_u_33']
            u12 = field_dict['_atom_site_aniso_u_12']
            u13 = field_dict['_atom_site_aniso_u_13']
            u23 = field_dict['_atom_site_aniso_u_23']
        except KeyError:
            return False
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        for value in values:
            try:
                atom = name_dict[self.index_to_str(label,value)]
            except KeyError:
                print(" CIF _atom_site_aniso_label, invalid key:", value)
                continue
            self.model.atom[atom].u_aniso = [
                self.index_to_float(u11,value),
                self.index_to_float(u22,value),
                self.index_to_float(u33,value),
                self.index_to_float(u12,value),
                self.index_to_float(u13,value),
                self.index_to_float(u23,value)
                ]
        return True

    def read_atom_site(self):
        for loop in self.loops:
            if self.read_atom_site_fract(None, loop.keys, loop.rows) or \
               self.read_atom_site_cartn(None, loop.keys, loop.rows):
                for loop in self.loops:
                    if self.read_atom_site_aniso(None, loop.keys, loop.rows):
                        break
                return True
        return False

    def read_struct_conn_(self, loop):
        try:
            type_id   = loop.get_col_idx('_struct_conn.conn_type_id')

            asym_id_1 = loop.get_col_idx('_struct_conn.ptnr1_auth_asym_id',
                                         '_struct_conn.ptnr1_label_asym_id')
            comp_id_1 = loop.get_col_idx('_struct_conn.ptnr1_auth_comp_id',
                                         '_struct_conn.ptnr1_label_comp_id')
            seq_id_1  = loop.get_col_idx('_struct_conn.ptnr1_auth_seq_id',
                                         '_struct_conn.ptnr1_label_seq_id')
            atom_id_1 = loop.get_col_idx('_struct_conn.ptnr1_label_atom_id')

            asym_id_2 = loop.get_col_idx('_struct_conn.ptnr2_auth_asym_id',
                                         '_struct_conn.ptnr2_label_asym_id')
            comp_id_2 = loop.get_col_idx('_struct_conn.ptnr2_auth_comp_id',
                                         '_struct_conn.ptnr2_label_comp_id')
            seq_id_2  = loop.get_col_idx('_struct_conn.ptnr2_auth_seq_id',
                                         '_struct_conn.ptnr2_label_seq_id')
            atom_id_2 = loop.get_col_idx('_struct_conn.ptnr2_label_atom_id')
        except KeyError:
            return False

        alt_id_1   = loop.get_col_idx_opt('_struct_conn.pdbx_ptnr1_label_alt_id')
        ins_code_1 = loop.get_col_idx_opt('_struct_conn.pdbx_ptnr1_pdb_ins_code')
        symm_1     = loop.get_col_idx_opt('_struct_conn.ptnr1_symmetry')
        alt_id_2   = loop.get_col_idx_opt('_struct_conn.pdbx_ptnr2_label_alt_id')
        ins_code_2 = loop.get_col_idx_opt('_struct_conn.pdbx_ptnr2_pdb_ins_code')
        symm_2     = loop.get_col_idx_opt('_struct_conn.ptnr2_symmetry')

        idxs_1 = [asym_id_1, comp_id_1, seq_id_1, ins_code_1, atom_id_1, alt_id_1]
        idxs_2 = [asym_id_2, comp_id_2, seq_id_2, ins_code_2, atom_id_2, alt_id_2]

        # atoms indexed by atomic identifiers
        atom_dict = dict(((a.chain, a.resn, a.resi,
            getattr(a, 'ins_code', ''), a.name, a.alt), i)
            for (i, a) in enumerate(self.model.atom))

        for row in loop.rows:
            if self.index_to_str(type_id, row).lower() != 'covale':
                # ignore non-covalent bonds (metalc, hydrog)
                continue
            if self.index_to_str(symm_1, row) != self.index_to_str(symm_2, row):
                # don't bond to symmetry mates
                continue
            key_1 = tuple(self.index_to_str(i, row) for i in idxs_1)
            key_2 = tuple(self.index_to_str(i, row) for i in idxs_2)
            try:
                index = [atom_dict[key_1], atom_dict[key_2]]
            except KeyError:
                print(" CIF _struct_conn, invalid keys:", row)
                continue
            bond = Bond()
            bond.index = index
            self.model.bond.append(bond)
        return True

    def read_struct_conn(self):
        '''
        Create bonds from STRUCT_CONN category
        '''
        for loop in self.loops:
            if self.read_struct_conn_(loop):
                return True
        return False

    def read_geom_bond_atom_site_labels(self,fields,field_dict,values):
        try:
            label_1 = field_dict.get('_geom_bond_atom_site_id_1', None)
            if label_1 is not None:
                label_2 = field_dict['_geom_bond_atom_site_id_2']
            else:
                label_1 = field_dict['_geom_bond_atom_site_label_1']
                label_2 = field_dict['_geom_bond_atom_site_label_2']
        except KeyError:
            return False

        symm_1 = field_dict.get('_geom_bond_site_symmetry_1', -1)
        symm_2 = field_dict.get('_geom_bond_site_symmetry_2', -1)

        # create index of atom name
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        for value in values:
            if self.index_to_str(symm_1, value) != self.index_to_str(symm_2, value):
                # don't bond to symmetry mates
                continue
            try:
                index = [name_dict[self.index_to_str(label_1, value)],
                         name_dict[self.index_to_str(label_2, value)]]
            except KeyError:
                print(" CIF _geom_bond_atom_site_label, invalid keys:", value)
                continue
            bond = Bond()
            bond.index = index
            bond.order = 1
            self.model.bond.append(bond)
        return True

    def read_geom_bond(self):
        '''
        Create bonds from GEOM_BOND category
        '''
        for loop in self.loops:
            if self.read_geom_bond_atom_site_labels(None, loop.keys, loop.rows):
                return True
        return False

    def read_chem_comp_bond_atom_ids(self,fields,field_dict,values):
        try:
            label_1 = field_dict['_chem_comp_bond_atom_id_1']
            label_2 = field_dict['_chem_comp_bond_atom_id_2']
        except KeyError:
            return False
        order_table = { 'sing' : 1, 'doub' : 2, 'trip' :3, 'delo': 4 }
        # create index of atom name
        cnt = 0
        name_dict = {}
        for atom in self.model.atom:
            if hasattr(atom,'name'):
                name_dict[atom.name] = cnt
            cnt = cnt + 1
        order = field_dict.get('_chem_comp_bond_value_order',None)
        for value in values:
            try:
                index = [name_dict[self.index_to_str(label_1, value)],
                         name_dict[self.index_to_str(label_2, value)]]
            except KeyError:
                print(" CIF _chem_comp_bond_atom_id, invalid keys:", value)
                continue
            bond = Bond()
            bond.index = index
            if order is not None:
                order_string = self.index_to_str(order,value).lower()
                bond.order = order_table.get(order_string[0:4],1)
            else:
                bond.order = 1
            self.model.bond.append(bond)
        # don't do distance-based bonding
        self.model.connect_mode = 1
        return True

    def read_chem_comp_bond(self):
        '''
        Create bonds from CHEM_COMP_BOND category.

        If successfull, don't do any distance-based bonding.
        '''
        for loop in self.loops:
            if self.read_chem_comp_bond_atom_ids(None, loop.keys, loop.rows):
                return True
        return False

    def read_ss_(self, loop, ss, ssrecords):
        '''
        Populate ssrecords dictionary with secondary structure records from
        STRUCT_CONF or STRUCT_SHEET_RANGE category.
        '''
        prefix = '_struct_conf' if ss == 'H' else '_struct_sheet_range'

        try:
            Beg_NDB_Strand_ID = loop.get_col_idx(prefix + '.beg_auth_asym_id',
                                                 prefix + '.beg_label_asym_id')
            Beg_NDB_res_num   = loop.get_col_idx(prefix + '.beg_auth_seq_id',
                                                 prefix + '.beg_label_seq_id')
            End_NDB_Strand_ID = loop.get_col_idx(prefix + '.end_auth_asym_id',
                                                 prefix + '.end_label_asym_id')
            End_NDB_res_num   = loop.get_col_idx(prefix + '.end_auth_seq_id',
                                                 prefix + '.end_label_seq_id')
        except KeyError:
            return False

        Beg_NDB_ins_code = loop.get_col_idx_opt(prefix + '.pdbx_beg_pdb_ins_code')
        End_NDB_ins_code = loop.get_col_idx_opt(prefix + '.pdbx_end_pdb_ins_code')

        idxs_1 = [Beg_NDB_Strand_ID, Beg_NDB_res_num, Beg_NDB_ins_code]
        idxs_2 = [End_NDB_Strand_ID, End_NDB_res_num, End_NDB_ins_code]

        for row in loop.rows:
            key = tuple(self.index_to_str(i, row) for i in idxs_1)
            ssrecords[key] = [ss, tuple(self.index_to_str(i, row) for i in idxs_2)]

        return True

    def read_ss(self):
        '''
        Read secondary structre (sheets and helices) from STRUCT_CONF and
        STRUCT_SHEET_RANGE categories and assign to CA atoms.
        '''
        ssrecords = {}

        for loop in self.loops:
            if self.read_ss_(loop, 'H', ssrecords):
                break

        for loop in self.loops:
            if self.read_ss_(loop, 'S', ssrecords):
                break

        if not ssrecords:
            return False

        atoms = self.model.atom

        for i, a in enumerate(atoms):
            if a.name != 'CA':
                continue

            try:
                ss, endkey = ssrecords[a.chain, a.resi,
                        getattr(a, 'ins_code', '')]
            except KeyError:
                continue

            for j in range(i, len(atoms)):
                aj = atoms[j]
                if aj.name != 'CA':
                    continue
                aj.ss = ss
                if (aj.chain, aj.resi, getattr(aj, 'ins_code', '')) == endkey:
                    break

        return True


class CIF:

    def __init__(self, fname, mode='r'):
        if mode not in ('r','pf'):
            print(" CIF: bad mode")
            return None
        if mode=='pf': # pseudofile
            contents = fname.read()
        else:
            try:
                from pymol.internal import file_read
                contents = file_read(fname)
            except ImportError:
                contents = open(fname, mode).read()
        self.datablocks_it = parse_cif(contents)

    def __iter__(self):
        return self

    def next(self):
        rec = CIFRec(next(self.datablocks_it))
        if rec.model.atom:
            return rec
        return next(self)

    __next__ = next

    def read(self):
        try:
            return next(self)
        except StopIteration:
            return None
