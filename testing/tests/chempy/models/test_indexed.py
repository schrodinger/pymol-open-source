from chempy import Atom, Bond
import pytest

from chempy.models import Indexed, Connected


@pytest.fixture
def indexed_model() -> Indexed:
    indexed_model = Indexed()

    for i in range(3):
        atom = Atom()
        atom.symbol = "H" if i == 0 else "He"
        atom.chain = "AB" if i == 0 else "BA"

        atom.coord = [float(i), 0.0, 0.0]
        indexed_model.add_atom(atom)

        bond = Bond()
        bond.index = [i, i]
        bond.order = i
        indexed_model.add_bond(bond)

    return indexed_model


def test_base_model_methods(indexed_model):
    assert indexed_model.get_residues() == [(0, 1), (1, 3)]
    assert indexed_model.get_coord_list() == [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
    ]
    assert indexed_model.get_mass() == 9.013144
    assert indexed_model.get_nuclear_charges() == 5

    assert indexed_model.index is None
    indexed_model.update_index()
    assert sorted(list(indexed_model.index.values())) == [0, 1, 2]


def test_get_implicit_mass(indexed_model):
    assert indexed_model.get_implicit_mass() == 10.021084


def test_list(indexed_model, capsys):
    indexed_model.list()
    captured = capsys.readouterr()
    assert captured.out == (
        "H  [0.0, 0.0, 0.0]\n"
        "He  [1.0, 0.0, 0.0]\n"
        "He  [2.0, 0.0, 0.0]\n"
        "[0, 0]\n"
        "[1, 1]\n"
        "[2, 2]\n"
    )


def test_get_min_max(indexed_model):
    assert indexed_model.get_min_max() == [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]


def test_indexed_merge(indexed_model):
    other_indexed_model = Indexed()

    atom = Atom()
    atom.symbol = "Ne"
    atom.chain = "F"

    atom.coord = [0.0, 0.0, 7.0]
    other_indexed_model.add_atom(atom)

    bond = Bond()
    bond.index = [1, 1]
    bond.order = 1
    other_indexed_model.add_bond(bond)

    indexed_model.update_index()
    assert sorted(list(indexed_model.index.values())) == [0, 1, 2]

    indexed_model.merge(other_indexed_model)

    assert sorted(list(indexed_model.index.values())) == [0, 1, 2, 3]
    assert len(indexed_model.atom) == 4
    assert len(indexed_model.bond) == 4

    assert other_indexed_model.index is None
    assert other_indexed_model.bond == []
    assert other_indexed_model.bond == []


def test_add_atom(indexed_model):
    new_atom = Atom()
    new_atom.symbol = "Ne"

    assert len(indexed_model.atom) == 3
    indexed_model.update_index()
    previos_last_atom = indexed_model.atom[-1]

    indexed_model.add_atom(new_atom)
    assert indexed_model.atom[-1] is new_atom
    assert indexed_model.atom[-2] is previos_last_atom

    assert len(indexed_model.atom) == 4


def test_delete_atom(indexed_model):

    assert len(indexed_model.atom) == 3
    assert len(indexed_model.bond) == 3
    assert indexed_model.atom[0].symbol == "H"

    indexed_model.delete_atom(0)

    assert len(indexed_model.atom) == 2
    assert len(indexed_model.bond) == 2
    assert indexed_model.atom[0].symbol == "He"
    assert indexed_model.bond[0].order == 1


def test_delete_list(indexed_model):

    assert len(indexed_model.atom) == 3
    assert len(indexed_model.bond) == 3
    assert indexed_model.atom[0].symbol == "H"

    indexed_model.delete_list([1, 2])

    assert len(indexed_model.atom) == 1
    assert len(indexed_model.bond) == 1
    assert indexed_model.atom[0].symbol == "H"
    assert indexed_model.bond[0].order == 0


def test_insert_atom(indexed_model):
    new_atom = Atom()
    new_atom.symbol = "Ne"

    assert len(indexed_model.atom) == 3
    previos_second_atom = indexed_model.atom[2]

    indexed_model.insert_atom(2, new_atom)

    assert indexed_model.atom[2] is new_atom
    assert previos_second_atom in indexed_model.atom
    assert indexed_model.atom.index(previos_second_atom) == 3
    assert indexed_model.index_atom(previos_second_atom) == 3
    assert len(indexed_model.atom) == 4


def test_index_atom(indexed_model):
    new_atom = Atom()
    new_atom.symbol = "Ne"

    indexed_model.add_atom(new_atom)

    assert indexed_model.index_atom(new_atom) == 3


def test_add_bond(indexed_model):
    assert len(indexed_model.atom) == 3
    assert len(indexed_model.bond) == 3
    assert indexed_model.bond[-1].order == 2

    new_bond = Bond()
    new_bond.index = [42, 42]
    new_bond.order = 42
    indexed_model.add_bond(new_bond)

    assert len(indexed_model.atom) == 3
    assert len(indexed_model.bond) == 4
    assert indexed_model.bond[-1].order == 42


def test_convert_to_connected(indexed_model):
    connected_model = indexed_model.convert_to_connected()
    assert isinstance(connected_model, Connected)
    assert len(connected_model.bond) == 3
    assert connected_model.bond[0] is connected_model.bond[0]
    # maybe a bug ?

    assert indexed_model.atom == []
    assert indexed_model.bond == []
    assert indexed_model.index == None


def test_get_internal_tuples(indexed_model):
    assert indexed_model.get_internal_tuples() == [(1,), (1, 1), (-1, 1, 1)]
