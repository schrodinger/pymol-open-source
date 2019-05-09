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

# mass calculation (even for implicit models)

implicit_valence = {
    'H'  :  {0:1,1:0},
    'C'  :  {0:4,1:3,2:2,3:1,4:0},
    'N'  :  {0:3,1:2,2:1,3:0},
    'O'  :  {0:2,1:1,2:0},
    'F'  :  {0:1,1:0},
    'Cl' :  {0:1,1:0},
    'CL' :  {0:1,1:0},
    'Br' :  {0:1,1:0},
    'BR' :  {0:1,1:0},
    'I'  :  {0:1,1:0},
    'S'  :  {0:2,1:2,2:0,3:1,4:0,5:1,6:0} # ambiguity?
    }

def implicit_mass(indexed):
    valence = [0]*len(indexed.atom)
    implicit = [0]*len(indexed.atom)

    for a in indexed.bond:
        ai0 = a.index[0]
        ai1 = a.index[1]
        valence[ai0] = valence[ai0] + 1
        valence[ai1] = valence[ai1] + 1
    c = 0
    for a in model.atom:
        valence[c] = valence[c] - a.formal_charge
        implicit[c] = implicit_valence[a.symbol][valence[c]]
    c = c + 1
