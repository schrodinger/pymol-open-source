'''
Experimental MMTF (Macromolecular Transmission Format) support
'''

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

try:
    from itertools import izip
    as_str = str
except ImportError:
    # python3
    izip = zip
    as_str = lambda s: s if isinstance(s, str) else s.decode()

#####################################################################

ss_map = {
    0: 'H', # pi helix
    1: 'L', # bend
    2: 'H', # alpha helix
    3: 'S', # extended
    4: 'H', # 3-10 helix
    5: 'S', # bridge
    6: 'L', # turn
    7: 'L', # coil
}

#####################################################################

def _to_chempy(data, use_auth=True):
    '''
    Construct a "chempy" model (molecule) from decoded MMTF data.
    '''
    from itertools import islice
    from chempy import models, Atom, Bond

    def add_bond(i1, i2, order, offset=0):
        bond = Bond()
        bond.order = order
        bond.index = [i1 + offset, i2 + offset]
        model.add_bond(bond)

    coord_iter = data.get_table_iter([
        'xCoordList',
        'yCoordList',
        'zCoordList',
    ])

    atom_iter = data.get_table_iter([
        'bFactorList',
        'occupancyList',
        'altLocList',
        'atomIdList',
    ], [0.0, 1.0, '', -1])

    group_iter = data.get_table_iter([
        'groupTypeList',
        'sequenceIndexList',
        'groupIdList',
        'insCodeList',
        'secStructList',
    ])

    chain_list_iter = enumerate(data.get_table_iter([
        'chainIdList',
        'chainNameList',
        'groupsPerChain',
    ]))

    groupList = data.get('groupList')

    symmetry = (
        data.get('unitCell', None),
        as_str(data.get('spaceGroup', '')),
    )

    model_output = []

    for n_chains in data.get_iter('chainsPerModel'):
        model = models.Indexed()
        model_output.append(model)

        if symmetry[0] is not None:
            model.cell, model.spacegroup = symmetry

        for (chain_idx, (segi, chain, n_groups)) in islice(chain_list_iter, n_chains):
            for (groupType, label_seq_id, auth_seq_id, ins_code, ss_info) in \
                    islice(group_iter, n_groups):

                group = groupList[groupType]
                resn = as_str(group[b'groupName'])

                group_bond_iter = izip(
                        group[b'bondAtomList'][0::2],
                        group[b'bondAtomList'][1::2],
                        group[b'bondOrderList'],
                        )

                offset = len(model.atom)
                for (i1, i2, order) in group_bond_iter:
                    add_bond(i1, i2, order, offset)

                group_atom_iter = izip(
                        group[b'atomNameList'],
                        group[b'elementList'],
                        group[b'formalChargeList'],
                        )

                for (name, elem, formal_charge) in group_atom_iter:
                    atom = Atom()

                    (atom.b, atom.q, atom.alt, atom.id) = next(atom_iter)

                    atom.coord = next(coord_iter)
                    atom.symbol = as_str(elem)
                    atom.name = as_str(name)
                    atom.resn = resn
                    atom.hetatm = label_seq_id == -1
                    atom.formal_charge = formal_charge
                    atom.segi = segi
                    atom.chain = chain
                    atom.ss = ss_map.get(ss_info, '')

                    if use_auth or label_seq_id is None:
                        atom.resi = auth_seq_id
                        atom.ins_code = ins_code or ''
                    else:
                        atom.resi = label_seq_id + 1

                    model.add_atom(atom)

    model_atom_max = 0
    model_atom_min = 0
    model_iter = iter(model_output)
    bondAtomList_iter = data.get_iter('bondAtomList')

    for order in data.get_iter('bondOrderList'):
        i1 = next(bondAtomList_iter)
        i2 = next(bondAtomList_iter)
        if i1 >= model_atom_max or i2 >= model_atom_max:
            model = next(model_iter)
            model_atom_min = model_atom_max
            model_atom_max += len(model.atom)
        add_bond(i1, i2, order, -model_atom_min)

    return model_output

#####################################################################

from .io import MmtfReader
MmtfReader.to_chempy = _to_chempy

#####################################################################
