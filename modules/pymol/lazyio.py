'''
Support for some less common file formats for PyMOL.

Copyright (c) Schrodinger, LLC.
'''

import os

import pymol
from pymol import cmd, CmdException

try:
    from lxml import etree
except ImportError:
    import xml.etree.ElementTree as etree


def load_pdbml(filename, object='', discrete=0, multiplex=1, zoom=-1,
        quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a PDBML formatted structure file
    '''
    from chempy import Atom, models
    from collections import defaultdict

    multiplex, discrete = int(multiplex), int(discrete)

    try:
        root = etree.fromstring(_self.file_read(filename))
        PDBxNS = root.tag.rstrip('datablock')
        atom_site_list = root.findall('.'
                '/' + PDBxNS + 'atom_siteCategory'
                '/' + PDBxNS + 'atom_site')
    except etree.XMLSyntaxError:
        raise CmdException("File doesn't look like XML")
    except etree.XPathEvalError:
        raise CmdException("XML file doesn't look like a PDBML file")

    if not atom_site_list:
        raise CmdException("no PDBx:atom_site nodes found in XML file")

    # state -> model dictionary
    model_dict = defaultdict(models.Indexed)

    # atoms
    for atom_site in atom_site_list:
        atom = Atom()
        atom.coord = [None, None, None]

        model_num = 1

        for child in atom_site:
            tag = child.tag

            if tag == PDBxNS + 'Cartn_x':
                atom.coord[0] = float(child.text)
            elif tag == PDBxNS + 'Cartn_y':
                atom.coord[1] = float(child.text)
            elif tag == PDBxNS + 'Cartn_z':
                atom.coord[2] = float(child.text)
            elif tag == PDBxNS + 'B_iso_or_equiv':
                atom.b = float(child.text)
            elif tag == PDBxNS + 'auth_asym_id':
                atom.chain = child.text or ''
            elif tag == PDBxNS + 'auth_atom_id':
                atom.name = child.text or ''
            elif tag == PDBxNS + 'auth_comp_id':
                atom.resn = child.text or ''
            elif tag == PDBxNS + 'auth_seq_id':
                atom.resi = child.text or ''
            elif tag == PDBxNS + 'label_alt_id':
                atom.resi = child.text or ''
            elif tag == PDBxNS + 'label_asym_id':
                atom.segi = child.text or ''
            elif tag == PDBxNS + 'label_atom_id':
                if not atom.name:
                    atom.name = child.text or ''
            elif tag == PDBxNS + 'label_comp_id':
                if not atom.resn:
                    atom.resn = child.text or ''
            elif tag == PDBxNS + 'label_seq_id':
                if not atom.resi:
                    atom.resi = child.text or ''
            elif tag == PDBxNS + 'label_entity_id':
                atom.custom = child.text or ''
            elif tag == PDBxNS + 'occupancy':
                atom.q = float(child.text)
            elif tag == PDBxNS + 'pdbx_PDB_model_num':
                model_num = int(child.text)
            elif tag == PDBxNS + 'type_symbol':
                atom.symbol = child.text or ''
            elif tag == PDBxNS + 'group_PDB':
                atom.hetatm = (child.text == 'HETATM')

        if None not in atom.coord:
            model_dict[model_num].add_atom(atom)

    # symmetry and cell
    try:
        node = root.findall('.'
                '/' + PDBxNS + 'cellCategory'
                '/' + PDBxNS + 'cell')[0]
        cell = [
            float(node.findall('./' + PDBxNS + a)[0].text)
            for a in [
                'length_a', 'length_b', 'length_c', 'angle_alpha',
                'angle_beta', 'angle_gamma'
            ]
        ]

        spacegroup = root.findall('.'
                '/' + PDBxNS + 'symmetryCategory'
                '/' + PDBxNS + 'symmetry'
                '/' + PDBxNS + 'space_group_name_H-M')[0].text
    except IndexError:
        cell = None
        spacegroup = ''

    # object name
    if not object:
        object = os.path.basename(filename).split('.', 1)[0]

    # only multiplex if more than one model/state
    multiplex = multiplex and len(model_dict) > 1

    # load models as objects or states
    for model_num in sorted(model_dict):
        if model_num < 1:
            print(" Error: model_num < 1 not supported")
            continue

        model = model_dict[model_num]
        model.connect_mode = 3

        if cell:
            model.cell = cell
            model.spacegroup = spacegroup

        if multiplex:
            oname = '%s_%04d' % (object, model_num)
            model_num = 1
        else:
            oname = object

        _self.load_model(model, oname,
                state=model_num, zoom=zoom, discrete=discrete)


def load_cml(filename, object='', discrete=0, multiplex=1, zoom=-1,
        quiet=1, _self=cmd):
    '''
DESCRIPTION

    Load a CML formatted structure file
    '''
    from chempy import Atom, Bond, models

    multiplex, discrete = int(multiplex), int(discrete)

    try:
        root = etree.fromstring(_self.file_read(filename))
    except etree.XMLSyntaxError:
        raise CmdException("File doesn't look like XML")

    if root.tag != 'cml':
        raise CmdException('not a CML file')

    molecule_list = root.findall('./molecule')

    if len(molecule_list) < 2:
        multiplex = 0
    elif not multiplex:
        discrete = 1

    for model_num, molecule_node in enumerate(molecule_list, 1):
        model = models.Indexed()

        atom_idx = {}

        for atom_node in molecule_node.findall('./atomArray/atom'):
            atom = Atom()
            atom.name = atom_node.get('id', '')

            if 'x3' in atom_node.attrib:
                atom.coord = [float(atom_node.get(a))
                        for a in ['x3', 'y3', 'z3']]
            elif 'x2' in atom_node.attrib:
                atom.coord = [float(atom_node.get(a))
                        for a in ['x2', 'y2']] + [0.0]
            else:
                print(' Warning: no coordinates for atom', atom.name)
                continue

            atom.symbol = atom_node.get('elementType', '')
            atom.formal_charge = int(atom_node.get('formalCharge', 0))
            atom_idx[atom.name] = len(model.atom)
            model.add_atom(atom)

        for bond_node in molecule_node.findall('./bondArray/bond'):
            refs = bond_node.get('atomsRefs2', '').split()
            if len(refs) == 2:
                bnd = Bond()
                bnd.index = [int(atom_idx[ref]) for ref in refs]
                bnd.order = int(bond_node.get('order', 1))
                model.add_bond(bnd)

        # object name
        if not object:
            object = os.path.basename(filename).split('.', 1)[0]

        # load models as objects or states
        if multiplex:
            oname = molecule_node.get('id') or _self.get_unused_name('unnamed')
            model_num = 1
        else:
            oname = object

        _self.load_model(model, oname,
                state=model_num, zoom=zoom, discrete=discrete)


def get_stlstr(binary=1, quiet=0, _self=cmd):
    '''
DESCRIPTION

    STL geometry export
    '''
    raise pymol.IncentiveOnlyException(
            "STL export not supported by this PyMOL build")


def read_stlstr(contents, object, state=0, zoom=-1, _self=cmd):
    '''
DESCRIPTION

    Load STL ASCII or binary content as a CGO object
    '''
    raise pymol.IncentiveOnlyException(
            "STL import not supported by this PyMOL build")


def read_collada(contents: bytes, oname, state=0, zoom=-1, *, _self=cmd):
    '''
DESCRIPTION

    Load COLLADA content as a CGO object
    '''
    raise pymol.IncentiveOnlyException(
            "COLLADA import not supported by this PyMOL build")


def get_mmtfstr(selection='all', state=1, _self=cmd):
    '''
DESCRIPTION

    DEPRECATED: Use cmd.get_bytes('mmtf')

    Export an atom selection to MMTF format.
    '''
    import simplemmtf

    try:
        # register PyMOL-specific spec extensions
        simplemmtf.levels[u'atom'][u'pymolReps'] = 0
        simplemmtf.levels[u'atom'][u'pymolColor'] = 0
        simplemmtf.encodingrules[u'pymolRepsList'] = (7, 0)
    except Exception as e:
        print(e)

    mmtfstr = simplemmtf.mmtfstr

    ss_map = {
        'H': 2,  # alpha helix
        'S': 3,  # extended
    }

    def callback(state, segi, chain, resv, resi, resn, name, elem,
            x, y, z, reps, color, alt, formal_charge, b, q, ss):
        atoms.append({
            u'modelIndex': state,
            u'chainId': mmtfstr(segi),
            u'chainName': mmtfstr(chain),
            u'groupId': resv,
            u'groupName': mmtfstr(resn),
            u'atomName': mmtfstr(name),
            u'element': mmtfstr(elem),
            u'coords': (x, y, z),
            u'altLoc': mmtfstr(alt),
            u'formalCharge': formal_charge,
            u'bFactor': b,
            u'occupancy': q,
            u'secStruct': ss_map.get(ss, -1),
            u'insCode': mmtfstr(resi.lstrip('0123456789')),
            u'pymolReps': reps,
            u'pymolColor': color,
        })

    atoms = []
    _self.iterate_state(state, selection,
            'callback(state, segi, chain, resv, resi, resn, name, elem, '
            'x, y, z, reps, color, alt, formal_charge, b, q, ss)',
            space={'callback': callback})

    bonds = _self.get_bonds(selection, state)

    d_out = simplemmtf.from_atoms(atoms, bonds)

    return d_out.encode()
