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

# Input/output facility

from chempy.pkl import PKL
from chempy.lst import LST
from chempy.pdb import PDB
from chempy.xyz import XYZ
from chempy.mol import MOL
from chempy.mol2 import MOL2
from chempy.arc import ARC
from chempy.gms import GMS
from chempy.mmd import MMD
from chempy.mae import MAE
from chempy.cc1 import CC1

pkl = PKL() # general object io
lst = LST() # general string-list io

# specific for indexed model objects:

pdb = PDB()
xyz = XYZ()
mol = MOL()
mol2 = MOL2()
arc = ARC()
gms = GMS() # OBSOLETE - PLEASE DO NOT USE
mmd = MMD()
mae = MAE()
cc1 = CC1()
