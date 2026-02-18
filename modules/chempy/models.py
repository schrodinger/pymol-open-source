# A* -------------------------------------------------------------------
# B* This file contains source code for the PyMOL computer program
# C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
# D* -------------------------------------------------------------------
# E* It is unlawful to modify or remove this copyright notice.
# F* -------------------------------------------------------------------
# G* Please see the accompanying LICENSE file for further information.
# H* -------------------------------------------------------------------
# I* Additional authors of this source file include:
# -*
# -*
# -*
# Z* -------------------------------------------------------------------

import abc
import copy
from typing import Generic, Optional, TypeVar, Union

import numpy as np

from chempy import Atom, Bond, Molecule, cpv, feedback

T = TypeVar("T", Bond, list[Bond])


class Base(abc.ABC, Generic[T]):

    def __init__(self):
        self.index: Optional[dict[int, int]] = None
        self.molecule = Molecule()
        self.atom: list[Atom] = []
        self.bond: list[T] = []

    @property
    @abc.abstractmethod
    def _bond_groups(self) -> list[list[Bond]]:
        raise NotImplementedError()

    def _handle_new_atom(self) -> None:
        """In case of Connected class we need to add an empty list to slef.bond"""
        pass

    def add_atom(self, atom: Atom) -> int:
        self.send_feedback(key="atoms", message=f": adding atom atom {atom.name}")
        self.atom.append(atom)

        self._handle_new_atom()

        index = len(self.atom) - 1
        if self.index:
            self.index[id(atom)] = index
        return index

    def delete_atom(self, index: int) -> None:
        self.send_feedback(key="atoms", message=f": deleting atom {index}.")

        # update index if it exists
        if self.index is not None:
            for key, value in self.index.items():
                if value > index:
                    self.index[key] = value - 1
            del self.index[id(self.atom[index])]

        # delete atom
        del self.atom[index]

        # delete bonds associated with this atom

        for bond_group in self._bond_groups:
            templist = [i for i, bond in enumerate(bond_group) if index in bond.index]
            for i, j in enumerate(templist):
                del bond_group[j - i]

        # re-index bond table
        for bond_group in self._bond_groups:
            for bond in bond_group:
                if bond.index[0] > index:
                    bond.index[0] -= 1
                if bond.index[1] > index:
                    bond.index[1] -= 1

    def insert_atom(self, index: int, atom: Atom) -> None:
        self.send_feedback(
            key="atoms", message=f": inserting atom {atom.name} before {index}."
        )
        self.atom.insert(index, atom)

        # re-index bond table
        for bond_pair in self._bond_groups:
            for bond in bond_pair:
                if bond.index[0] >= index:
                    bond.index[0] = bond.index[0] + 1
                if bond.index[1] >= index:
                    bond.index[1] = bond.index[1] + 1

        # update index if it exists
        if self.index is not None:

            for key, value in self.index.items():
                if value >= index:
                    self.index[key] = value + 1

            self.index[id(atom)] = index

    def reset(self):
        self.index = None
        self.molecule = Molecule()
        self.atom.clear()
        self.bond.clear()

    def update_index(self) -> None:
        self.send_feedback(key="verbose", message="updating indexes...")
        self.index = {id(atom): index for index, atom in enumerate(self.atom)}

    def get_residues(self) -> list[tuple[int, int]]:
        result = []

        if not self.atom:
            return result

        last = self.atom[0]
        counter = 0
        start = 0

        for atom in self.atom:
            if not atom.in_same_residue(last):
                result.append((start, counter))
                start = counter
                last = atom
            counter += 1
        if counter - start > 1:
            result.append((start, counter))

        return result

    def get_coord_list(self) -> list[tuple[float, float, float]]:
        return [atom.coord for atom in self.atom]

    def get_mass(self) -> float:
        return sum(atom.get_mass() for atom in self.atom)

    def get_nuclear_charges(self):
        """Return the sum of nuclear charges of all atoms in a molecule."""
        return sum(atom.get_number() for atom in self.atom)

    def assign_names(self, preserve: bool = True) -> None:
        names = set()
        counter = {}

        for atom in self.atom:
            if preserve:
                if atom.has("name"):
                    names.add(atom.name)
            else:
                if hasattr(atom, "name"):
                    del atom.name

        for atom in self.atom:

            if atom.has("name"):
                continue

            c = counter[atom.symbol] if atom.symbol in counter else 1
            name = f"{atom.symbol}{c}"

            while name in names:
                name = f"{atom.symbol}{c}"
                c += 1

            counter[atom.symbol] = c
            atom.name = name

            names.add(name)

    def send_feedback(self, key: str, message: str) -> None:
        if feedback[key]:
            print(f" {self.__class__}): {message}")


class Indexed(Base[Bond]):

    @property
    def _bond_groups(self) -> list[list[Bond]]:
        return [self.bond]

    def get_min_max(self) -> Union[list[list[float]], list[tuple[float, float, float]]]:

        if not self.atom:
            return [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

        coords = np.array([atom.coord for atom in self.atom])

        min_coord = np.min(coords, axis=0).tolist()
        max_coord = np.max(coords, axis=0).tolist()

        return [min_coord, max_coord]

    def merge(self, other: "Indexed") -> None:
        """steals atom objects from 'other' and resets 'other'"""
        self.send_feedback(key="actions", message="merging models...")

        self.atom.extend(other.atom)

        for bond in other.bond:
            bond.index[0] += len(self.atom)
            bond.index[1] += len(self.atom)
            self.bond.append(bond)

        other.reset()

        if self.index:
            self.update_index()

    def delete_list(self, indexes: list[int]) -> None:
        """delete a list of indexed atoms."""

        if not indexes:
            return

        self.send_feedback(key="atoms", message=f": deleting atoms {indexes}")

        indexes_copy = copy.deepcopy(indexes)
        indexes_copy.sort()
        indexes_copy.reverse()

        # generate cross-reference tables
        shft = 0
        old_to_new = {}
        next = indexes_copy.pop()

        for i in range(len(self.atom)):
            if i == next:
                old_to_new[i] = -1
                next = indexes_copy.pop() if indexes_copy else None
                shft = shft - 1
            else:
                old_to_new[i] = i + shft

        if shft == 0:
            return

        # delete atoms
        new_atoms = []
        for i in range(len(self.atom)):
            if old_to_new[i] >= 0:
                new_atoms.append(self.atom[i])
        self.atom = new_atoms

        # delete bonds
        new_bonds: list[Bond] = []
        for bond in self.bond:
            b0 = bond.index[0]
            b1 = bond.index[1]
            if old_to_new[b0] >= 0 and old_to_new[b1] >= 0:
                bond.index[0] = old_to_new[b0]
                bond.index[1] = old_to_new[b1]
                new_bonds.append(bond)
        self.bond = new_bonds

        # update index if it exists
        if self.index is not None:
            self.index = {self.index[id(atom)]: i for i, atom in enumerate(self.atom)}

    def index_atom(self, atom: Atom):
        return self.atom.index(atom) if atom in self.atom else -1

    def add_bond(self, bond: Bond) -> None:
        self.send_feedback(
            key="bonds", message=f": adding bond ({bond.index[0]}, {bond.index[1]})"
        )
        self.bond.append(bond)

    def remove_bond(self, index: int) -> None:
        self.send_feedback(key="bonds", message=f": removing bond {index}")
        del self.bond[index]

    def convert_to_connected(self) -> "Connected":
        self.send_feedback(key="verbose", message=": converting to connected model...")

        model = Connected()
        model.molecule = self.molecule
        model.atom = self.atom
        model.bond = []
        model.index = None

        for _ in model.atom:
            model.bond.append([])

        for bond in self.bond:
            model.bond[bond.index[0]].append(bond)  # note two refs to same object
            model.bond[bond.index[1]].append(bond)  # note two refs to same object

        self.reset()

        return model

    def from_molobj(self, molobj: Molecule) -> None:
        self.reset()

        if len(molobj.title):
            self.molecule.title = molobj.title

        if len(molobj.comments):
            self.molecule.comments = molobj.comments

        self.molecule.chiral = molobj.chiral
        self.molecule.dim_code = molobj.dim_code

        for a in molobj.atom:
            at = Atom()
            at.symbol = a.symbol
            at.name = a.name

            if a.resn != Atom.defaults["resn"]:
                at.resn = a.resn

            if a.resn_code != Atom.defaults["resn_code"]:
                at.resn_code = a.resn_code

            at.resi = a.resi
            at.resi_number = a.resi_number
            at.b = a.b
            at.q = a.q
            at.alt = a.alt
            at.hetatm = a.hetatm

            if a.segi != Atom.defaults["segi"]:
                at.segi = a.segi

            if a.chain != Atom.defaults["chain"]:
                at.chain = a.chain

            at.color_code = a.color_code
            at.coord = a.coord
            at.formal_charge = a.formal_charge
            at.partial_charge = a.partial_charge

            if a.numeric_type != -99:
                at.numeric_type = a.numeric_type

            if a.text_type != "UNKNOWN":
                at.text_type = a.text_type

            at.stereo = a.stereo

            if hasattr(a, "flags"):
                at.flags = a.flags
            if hasattr(a, "vdw"):
                at.vdw = a.vdw

            self.atom.append(at)

        for b in molobj.bond:
            bnd = Bond()
            bnd.index = [b.atom[0], b.atom[1]]
            bnd.order = b.order
            bnd.stereo = b.stereo
            self.bond.append(bnd)

    def sort(self):
        self.send_feedback(key="verbose", message=f": sorting...")

        if self.index is None:
            self.update_index()

        old_index = self.index
        self.atom.sort()
        self.update_index()

        xref = {}
        new_index = self.index

        for a in list(new_index.keys()):
            xref[old_index[a]] = new_index[a]

        for b in self.bond:
            b.index[0] = xref[b.index[0]]
            b.index[1] = xref[b.index[1]]

    def get_internal_tuples(self) -> list:
        """
        generates raw atom sets needed to construct an internal coordinate
        description of the molecule
        get a connected version too
        """
        model_copy = copy.deepcopy(self).convert_to_connected()
        center = [0.0, 0.0, 0.0]
        to_go = len(self.atom)
        done = {}

        if to_go < 3:
            return [(0), (1, 0)]

        # get center of molecule
        for atom in self.atom:
            center = cpv.add(center, atom.coord)
        center = cpv.scale(center, 1.0 / len(self.atom))

        # find most central multivalent atom
        min_a = -1
        for counter, atom in enumerate(self.atom):
            if len(model_copy.bond[counter]) > 1:  # must have at least two neighbors
                d = cpv.distance(atom.coord, center)
                if min_a < 0:
                    min_d = d
                    min_a = counter
                elif d < min_d:
                    #    ^ not exists
                    min_d = d
                    min_a = counter

        # make this our first atom
        fst = min_a
        z_set = [(fst,)]
        done[fst] = 1
        to_go = to_go - 1

        # for the second atom, try to choose different multivalent neighbor
        nxt = -1
        for bond in model_copy.bond[fst]:
            neighbor = bond.index[0]
            if neighbor == fst:
                neighbor = bond.index[1]
            if len(model_copy.bond[neighbor]) > 1:
                nxt = neighbor
                break

        # safety, choose any neighbor
        if nxt < 0:
            neighbor = bond.index[0]
            if neighbor == fst:
                neighbor = bond.index[1]
            nxt = neighbor
        z_set.append((nxt, fst))
        done[nxt] = 1
        to_go = to_go - 1

        # for the third atom, choose a different multivalent neighbor
        trd = -1
        for bond in model_copy.bond[fst]:
            neighbor = bond.index[0]
            if neighbor == fst:
                neighbor = bond.index[1]
            if len(model_copy.bond[neighbor]) > 1:
                if neighbor not in done:
                    trd = neighbor
                    break

        # safety, choose any unchosen neighbor
        if trd < 0:
            for bond in model_copy.bond[fst]:
                neighbor = bond.index[0]
                if neighbor == fst:
                    neighbor = bond.index[1]
                if neighbor not in done:
                    trd = neighbor
                    break

        z_set.append((trd, fst, nxt))
        done[trd] = 1
        to_go = to_go - 1
        if to_go:
            # now find all torsions in the molecule
            tors = {}
            for bond in self.bond:  # use bond as center of torsion
                a1 = bond.index[0]
                a2 = bond.index[1]
                for counter in model_copy.bond[a1]:
                    a0 = counter.index[0]
                    if a0 not in (a1, a2):  # outside atom
                        for d in model_copy.bond[a2]:
                            a3 = d.index[0]
                            if a3 not in (a0, a1, a2):  # outside atom
                                if a0 < a3:
                                    to = (a0, a1, a2, a3)
                                else:
                                    to = (a3, a2, a1, a0)
                                tors[to] = 1
                            a3 = d.index[1]
                            if a3 not in (a0, a1, a2):  # outside atom
                                if a0 < a3:
                                    to = (a0, a1, a2, a3)
                                else:
                                    to = (a3, a2, a1, a0)
                                tors[to] = 1
                    a0 = counter.index[1]
                    if a0 not in (a1, a2):  # outside atom
                        for d in model_copy.bond[a2]:
                            a3 = d.index[0]
                            if a3 not in (a0, a1, a2):  # outside atom
                                if a0 < a3:
                                    to = (a0, a1, a2, a3)
                                else:
                                    to = (a3, a2, a1, a0)
                                tors[to] = 1
                            a3 = d.index[1]
                            if a3 not in (a0, a1, a2):  # outside atom
                                if a0 < a3:
                                    to = (a0, a1, a2, a3)
                                else:
                                    to = (a3, a2, a1, a0)
                                tors[to] = 1
            if len(tors):
                # choose remaining atoms based on existing atoms using torsion
                while to_go:
                    for tor in list(tors.keys()):
                        a0 = tor[0]
                        a1 = tor[1]
                        a2 = tor[2]
                        a3 = tor[3]
                        dh0 = a0 in done
                        dh1 = a1 in done
                        dh2 = a2 in done
                        dh3 = a3 in done
                        if (not dh0) and dh1 and dh2 and dh3:
                            z_set.append((a0, a1, a2, a3))
                            done[a0] = 1
                            to_go = to_go - 1
                        elif dh0 and dh1 and dh2 and (not dh3):
                            z_set.append((a3, a2, a1, a0))
                            done[a3] = 1
                            to_go = to_go - 1
            else:  # for molecules with no torsions (dichloromethane, etc.)
                # we have to generate torsions which include one virtual
                # bond
                for bond in model_copy.bond[fst]:
                    neighbor = bond.index[0]
                    if neighbor in [fst, nxt, trd]:
                        neighbor = bond.index[1]
                    if neighbor not in done:
                        z_set.append((neighbor, trd, fst, nxt))
                        to_go = to_go - 1
                        done[neighbor] = 1
        return z_set

    def list(self):
        for atom in self.atom:
            print(atom.symbol, atom.name, atom.coord)
        for bond in self.bond:
            print(bond.index)

    def get_implicit_mass(self):
        """mass calculation for implicit models"""

        valence = [0] * len(self.atom)

        for bond in self.bond:
            v = 1.5 if bond.order == 4 else bond.order
            valence[bond.index[0]] += v
            valence[bond.index[1]] += v

        h_count = sum(atom.get_free_valence(v) for (atom, v) in zip(self.atom, valence))

        hydrogen = Atom()
        hydrogen.symbol = "H"
        return self.get_mass() + hydrogen.get_mass() * h_count


class Connected(Base[list[Bond]]):

    @property
    def _bond_groups(self) -> list[list[Bond]]:
        return self.bond

    def _handle_new_atom(self) -> None:
        self.bond.append([])

    def convert_to_indexed(self) -> Indexed:
        self.send_feedback(key="actions", message=": converting to indexed model...")

        indexed = Indexed()
        indexed.atom = self.atom
        indexed.molecule = self.molecule
        for counter, bond_pair in enumerate(self.bond):
            for bond in bond_pair:
                if bond.index[0] == counter:
                    indexed.bond.append(bond)
        # TODO: remove next line, otherwise indexed model will be empty
        # self.reset()
        return indexed

    def sort(self):
        self.send_feedback(key="verbose", message=f": sorting...")

        if self.index is None:
            self.update_index()

        old_index = self.index
        self.atom.sort()
        self.update_index()

        xref = {}
        new_index = self.index

        for a in new_index.keys():
            xref[old_index[a]] = new_index[a]

        new_bond = [None] * len(self.atom)
        c = 0
        tmp_list = []
        for bond_pair in self.bond:
            for bond in bond_pair:
                if c == bond.index[0]:
                    tmp_list.append(bond)
            new_bond[xref[c]] = bond_pair
            c = c + 1

        for b in tmp_list:
            b.index[0] = xref[b.index[0]]
            b.index[1] = xref[b.index[1]]

        self.bond = new_bond
