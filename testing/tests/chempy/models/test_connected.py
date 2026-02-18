from chempy import Atom, Bond
import pytest

from chempy.models import Indexed, Connected


@pytest.fixture
def connected_model() -> Connected:
    connected_model = Connected()

    for i in range(3):
        atom = Atom()
        atom.symbol = "H" if i == 0 else "He"
        atom.chain = "AB" if i == 0 else "BA"

        atom.coord = [float(i), 0.0, 0.0]
        connected_model.add_atom(atom)

    for i in range(len(connected_model.atom)):
        bond = Bond()
        bond.index = [i, i]
        bond.order = i
        # note two refs to same object
        connected_model.bond[bond.index[0]].append(bond)
        connected_model.bond[bond.index[1]].append(bond)
    return connected_model


def test_base_model_methods(connected_model):
    assert len(connected_model.atom) == 3
    assert len(connected_model.bond) == 3
    assert connected_model.get_residues() == [(0, 1), (1, 3)]
    assert connected_model.get_coord_list() == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ]
    assert connected_model.get_mass() == 9.013144
    assert connected_model.get_nuclear_charges() == 5

    assert connected_model.index is None
    connected_model.update_index()
    assert sorted(list(connected_model.index.values())) == [0, 1, 2]


def test_add_atom(connected_model):
    new_atom = Atom()
    new_atom.symbol = "Ne"

    assert len(connected_model.atom) == 3
    connected_model.update_index()
    previos_last_atom = connected_model.atom[0]

    connected_model.add_atom(new_atom)
    assert connected_model.atom[-1] is new_atom
    assert connected_model.atom[-2] is not previos_last_atom

    assert len(connected_model.atom) == 4


def test_delete_atom(connected_model):

    assert len(connected_model.atom) == 3
    assert len(connected_model.bond) == 3
    assert connected_model.atom[0].symbol == "H"

    connected_model.delete_atom(0)

    assert len(connected_model.atom) == 2
    assert len(connected_model.bond) == 3
    assert connected_model.bond[0] == []
    assert connected_model.atom[0].symbol == "He"
    assert connected_model.bond[1][0].order == 1


def test_insert_atom(connected_model):
    new_atom = Atom()
    new_atom.symbol = "Ne"

    assert len(connected_model.atom) == 3
    previos_second_atom = connected_model.atom[2]

    connected_model.insert_atom(2, new_atom)

    assert connected_model.atom[2] is new_atom
    assert previos_second_atom in connected_model.atom
    assert connected_model.atom.index(previos_second_atom) == 3
    assert len(connected_model.atom) == 4


def test_convert_to_indexed(connected_model):
    indexed_model = connected_model.convert_to_indexed()
    assert isinstance(indexed_model, Indexed)
    assert len(indexed_model.atom) == 3
    # maybe a bug ?
    assert len(indexed_model.bond) == 6
