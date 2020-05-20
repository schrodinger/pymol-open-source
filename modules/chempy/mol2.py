#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright Schrodinger LLC.
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

from chempy import Storage

# see layer2/AtomInfo.h
cAtomFlag_polymer      = 0x08000000
cAtomFlag_solvent      = 0x10000000

class MOL2(Storage):

    _bondTypes = { 1         : "1",
                   2         : "2",
                   3         : "3",
                   4       : "ar",
                   0                : "nc",
                 }

    def toList(self,model,**kwargs):
        buf = ["# created with PyMOL\n"]

        # RTI MOLECULE
        buf.append("@<TRIPOS>MOLECULE\n")
        buf.append(model.molecule.title + "\n")

        # defer until number of substructures known
        buf_i_counts = len(buf)
        buf.append(None)

        buf.append("SMALL\n")
        buf.append("USER_CHARGES\n")

        no_text_type_count = 0
        at_prev = None
        residues = []

        # RTI ATOM
        # atom_id atom_name x y z atom_type [subst_id
        #   [subst_name [charge [status_bit]]]]
        buf.append("@<TRIPOS>ATOM\n")
        for atom_id, at in enumerate(model.atom, 1):
            resn = at.resn or "UNK"
            subst_name = resn + at.resi

            if not (at_prev and at_prev.in_same_residue(at)):
                residues.append([subst_name, atom_id,
                    at.chain or at.segi or '****',
                    resn, at.flags])
                at_prev = at

            text_type = at.text_type

            if not text_type:
                no_text_type_count += 1
                text_type = at.symbol or "Any"

            buf.append("%d\t%4s\t%.3f\t%.3f\t%.3f\t%2s\t%d\t%s\t%.3f\t%s\n" %
                           (atom_id,
                               at.name or at.symbol or "X",
                               at.coord[0],at.coord[1],at.coord[2],
                            text_type,
                            len(residues), subst_name,
                            at.partial_charge,
                            'WATER' if (at.flags & cAtomFlag_solvent) else '',
                            ))

        if no_text_type_count > 0:
            print(" Warning: %d atoms missing 'text_type', using element symbol instead."
                " Hint: use cmd.assign_atom_types() to assing MOL2 atom types.")

        # RTI BOND
        # bond_id origin_atom_id target_atom_id bond_type [status_bits]
        buf.append("@<TRIPOS>BOND\n")
        for b, bo in enumerate(model.bond):
            bOrder = MOL2._bondTypes[bo.order]
            buf.append("%d %d %d %s\n" % (b,
              1 + bo.index[0],
              1 + bo.index[1], bOrder))

        # RTI SUBSTRUCTURE
        # subst_id subst_name root_atom [subst_type [dict_type
        #   [chain [sub_type [inter_bonds [status [comment]]]]]]]
        buf.append("@<TRIPOS>SUBSTRUCTURE\n")
        for subst_id, res in enumerate(residues, 1):
            buf.append('%d\t%s\t%d\t%s\t1 %s\t%s\n' % (subst_id,
                res[0], res[1],
                'RESIDUE' if (res[4] & cAtomFlag_polymer) else 'GROUP',
                res[2], res[3]))

        buf[buf_i_counts] = "%d %d %d\n" % (model.nAtom, model.nBond, len(residues))

        return buf
