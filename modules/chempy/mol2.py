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

from __future__ import print_function

from chempy import Storage

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
        buf.append("%d %d\n" % (model.nAtom, model.nBond))
        buf.append("SMALL\n")
        buf.append("USER_CHARGES\n")

        no_text_type_count = 0

        # RTI ATOM
        buf.append("@<TRIPOS>ATOM\n")
        for at in model.atom:
            text_type = at.text_type

            if not text_type:
                no_text_type_count += 1
                text_type = at.symbol or "Any"

            buf.append("%d\t%4s\t%.3f\t%.3f\t%.3f\t%2s\t%d\t%s%s\t%.3f\n" %
                           (at.index,
                               at.name or at.symbol or "X",
                               at.coord[0],at.coord[1],at.coord[2],
                            text_type,
                            at.resi_number,
                            at.resn or "UNK",
                            at.resi,
                            at.q))

        if no_text_type_count > 0:
            print(" Warning: %d atoms missing 'text_type', using element symbol instead."
                " Hint: use cmd.assign_atom_types() to assing MOL2 atom types.")

        # RTI BOND
        buf.append("@<TRIPOS>BOND\n")
        for b, bo in enumerate(model.bond):
            bOrder = MOL2._bondTypes[bo.order]
            buf.append("%d %d %d %s\n" % (b,
              1 + bo.index[0],
              1 + bo.index[1], bOrder))

        return buf
