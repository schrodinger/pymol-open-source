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

from chempy.models import Indexed
from chempy import Storage,Atom,Bond,feedback

import string

class MOL2(Storage):

    def __init__(self,**kwargs):
        if 'cmd' in kwargs:
            self.cmd = kwargs['cmd']
        else:
            self.cmd = None

    # RTIs
    _fields = { "alt_type"               : "@<TRIPOS>ALT_TYPE\n", 
                "anchor_atom"            : "@<TRIPOS>ANCHOR_ATOM\n", 
                "associated_annotation"  : "@<TRIPOS>ASSOCIATED_ANNOTATION\n", 
                "atom"                   : "@<TRIPOS>ATOM\n", 
                "bond"                   : "@<TRIPOS>BOND\n", 
                "center_of_mass"         : "@<TRIPOS>CENTER_OF_MASS\n", 
                "centroid"               : "@<TRIPOS>CENTROID\n", 
                "comment"                : "@<TRIPOS>COMMENT\n", 
                "crysin"                 : "@<TRIPOS>CRYSIN\n", 
                "dict"                   : "@<TRIPOS>DICT\n", 
                "data_file"              : "@<TRIPOS>DATA_FILE\n", 
                "extension_point"        : "@<TRIPOS>EXTENSION_POINT\n", 
                "ff_pbc"                 : "@<TRIPOS>FF_PBC\n", 
                "ffcon_angle"            : "@<TRIPOS>FFCON_ANGLE\n", 
                "ffcon_dist"             : "@<TRIPOS>FFCON_DIST\n", 
                "ffcon_multi"            : "@<TRIPOS>FFCON_MULTI\n", 
                "ffcon_range"            : "@<TRIPOS>FFCON_RANGE\n", 
                "ffcon_torsion"          : "@<TRIPOS>FFCON_TORSION\n", 
                "line"                   : "@<TRIPOS>LINE\n", 
                "lsplane"                : "@<TRIPOS>LSPLANE\n", 
                "molecule"               : "@<TRIPOS>MOLECULE\n", 
                "normal"                 : "@<TRIPOS>NORMAL\n", 
                "qsar_align_rule"        : "@<TRIPOS>QSAR_ALIGN_RULE\n", 
                "ring_closure"           : "@<TRIPOS>RING_CLOSURE\n", 
                "rotatable_bond"         : "@<TRIPOS>ROTATABLE_BOND\n", 
                "search_dist"            : "@<TRIPOS>SEARCH_DIST\n", 
                "search_options"         : "@<TRIPOS>SEARCH_OPTIONS\n", 
                "set"                    : "@<TRIPOS>SET\n", 
                "substructure"           : "@<TRIPOS>SUBSTRUCTURE\n", 
                "u_feat"                 : "@<TRIPOS>U_FEAT\n", 
                "unity_atom_attr"        : "@<TRIPOS>UNITY_ATOM_ATTR\n", 
                "unity_bond_attr"        : "@<TRIPOS>UNITY_BOND_ATTR\n" }

    _molType = { "small"            : "SMALL\n",
                 "bio"              : "BIOPOLYMER\n",
                 "prot"             : "PROTEIN\n",
                 "nuc"              : "NUCLEIC_ACID\n",
                 "sacc"             : "SACCHARIDE\n" }

    _chargeType = { "none"          : "NO_CHARGES\n",
                    "del"           : "DEL_RE\n",
                    "gast"          : "GASTEIGER\n",
                    "gast_h"        : "GAST_HUCK\n",
                    "huck"          : "HUCKEL\n",
                    "pull"          : "PULLMAN\n",
                    "gauss80"       : "GAUSS80_CHARGES\n",
                    "ampac"         : "AMPAC_CHARGES\n",
                    "mull"          : "MULLIKEN_CHARGES\n",
                    "dict"          : "DICT_ CHARGES\n",
                    "mmff94"        : "MMFF94_CHARGES\n",
                    "user"          : "USER_CHARGES\n" }

    _bondTypes = { 1         : "1",
                   2         : "2",
                   3         : "3",
                   "amide"          : "am",
                   4       : "ar",
                   "dummy"          : "du",
                   "unknown"        : "un",
                   "not_connected"  : "nc" }
                   

    def fromList(self,molList):
        pass
    
    def toList(self,model,**kwargs):
        molList = []

        if 'state' in kwargs:
            state = kwargs['state']
        else:
            state = None
        if 'selection' in kwargs:
            sel = kwargs['selection']
        else:
            sel = None

        f = MOL2._fields
        c = MOL2._chargeType
        m = MOL2._molType
        n = "\n"
        
        molList.append("# created with PyMOL")
        if model.molecule.comments!='':
            molList.append("# COMMENTS:"+n)
            molList.append("# " + model.molecule.comments)

        # RTI MOLECULE
        molList.append(n+f["molecule"])
        molList.append(model.molecule.title+n)
        subst=feat=sets=0
        molList.append("%d\t%d\t%d\t%d\t%d\n" % (model.nAtom,model.nBond,subst,feat,sets))
        # TODO: Guess this from the user's selection
        mKey="prot"
        ## if self.cmd!=None and state!=None and sel!=None:
        ##     nPoly = self.cmd.count_atoms("poly and (%s and state %s)" % (sel,state))
        ##     nOrg = self.cmd.count_atoms("org and (%s and state %s)" % (sel,state))
        ##     nIno = self.cmd.count_atoms("inorganic and (%s and state %s)" % (sel,state))
        ##     nNuc = self.cmd.count_atoms("(resn DA+DG+DC+DT+A+C+G+U) and (%s and state %s)" (sel,state))
        ##     # not too happy w/this
        ##     # - poly
        ##     if nPoly==0:
        ##         # - nuc
        ##         if nNuc==0:
        ##             # - org
        ##             if nOrg==0:
        ##                 # - ino
        ##                 if nIno==0:
        ##                     # what are you?!
        ##                     # you have no polymer, no organic, no nucleic, no inorganic
        ##                     # default, is wrong, but so are the others
        ##                     mKey = "prot"
        ##                 # ino - all
        ##                 else:
        ##                     mKey = "small"
        ##             # org
        ##             else:
        ##                 mKey = "small"
        ##         # nuc
        ##         else:
        ##             if nIno==0 and nOrg==0:
        ##                 mKey = "nuc"
        ##             else:
        ##                 mKey = "bio"
        ##     # poly
        ##     else:
        ##         mKey="prot"
            
        molList.append(m[mKey])
        # TODO: Guess this
        molList.append(c["user"])

        # RTI ATOM
        molList.append(f["atom"])
        for a in range(len(model.atom)):
            at = model.atom[a]
            molList.append("%d\t%4s\t%.3f\t%.3f\t%.3f\t%2s\t%.3f\n" %
                           (at.index,at.name,at.coord[0],at.coord[1],at.coord[2],
                            at.text_type, at.q))

        # RTI BOND
        molList.append(f["bond"])
        for b in range(len(model.bond)):
            bo = model.bond[b]
            bOrder = MOL2._bondTypes[bo.order]
            molList.append("%d %d %d %s\n" % (b,1+bo.index[0],1+bo.index[1],str(bOrder)))
        molList.append("\n")
        return molList

    def strToFile(self,dat,fname,**params):
        if feedback['io']:
            print ' chempy: writing mol2 to file "%s".' % fname
        fp = open(fname,'w')
        result = fp.writelines(dat)
        fp.close()
        
