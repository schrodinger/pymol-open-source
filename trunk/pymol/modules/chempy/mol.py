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

from chempy.models import Indexed
from chempy import Storage,Atom,Bond

import string

class MOL(Storage):

    def fromList(self,molList):

        model = Indexed()

        # read header information
        model.molecule.title = string.strip(molList[0])
        model.molecule.dim_code = string.strip(molList[1][20:22])
        model.molecule.comments = string.strip(molList[2])
        try:
            model.molecule.chiral = int(molList[3][12:15])
        except:
            model.molecule.chiral = 0
        nAtom = int(molList[3][0:3])
        nBond = int(molList[3][3:6])

        # read atoms
        nameDict = {}
        irec = 4
        cnt = 0
        for a in range(nAtom):
            at = Atom()
            at.index = cnt
            at.coord = [float(molList[irec][0:10]), 
                float(molList[irec][10:20]),float(molList[irec][20:30])]
            at.symbol = string.strip(molList[irec][31:33])
            try:
                at.stereo = int(molList[irec][39:42])
            except:
                at.stereo = 0
            chg=int(molList[irec][36:39])
            if chg>0: chg=4-chg
            at.formal_charge = chg
            model.atom.append(at)
            irec = irec + 1
            cnt = cnt + 1

            # read bonds
        for a in range(nBond):
            bnd = Bond()
            bnd.index = [ int(molList[irec][0:3])-1,int(molList[irec][3:6])-1 ]
            bnd.order = int(molList[irec][6:9])
            try:
                bnd.stereo = int(molList[irec][9:12])
            except:
                bnd.stereo = 0
            model.bond.append(bnd)
            irec = irec+1

            # obtain formal charges from M  CHG record
        while molList[irec][0:6]!='M  END':
            if molList[irec][0:6]=='M  CHG':
                cl = string.split(string.strip(molList[irec][6:]))
                cll = int(cl[0])*2
                a=1
                while a<=cll:
                    model.atom[int(cl[a])-1].formal_charge=int(cl[a+1])
                    a=a+2
            irec =irec+1
            if irec >= len(molList): break

        return model

#------------------------------------------------------------------------------
    def toList(self,model):

        molList = []

        # write header records
        molList.append(model.molecule.title+"\n")
        molList.append("  ChemPy            %2s                             0\n" %
                 model.molecule.dim_code)
        molList.append(model.molecule.comments+"\n")
        molList.append("%3d%3d  0  0  %1d  0  0  0  0  0999 V2000\n" %
                            (model.nAtom, model.nBond, model.molecule.chiral))

        # write atom records
        for a in model.atom:
            chg = a.formal_charge
            if chg!=0: chg=4-chg
            molList.append("%10.4f%10.4f%10.4f %-3s 0  %1d  %1d  0  0  0  0  0  0  0  0  0\n" % \
                            (a.coord[0], a.coord[1], a.coord[2], a.symbol, chg, a.stereo))

            # write bond records
        for b in model.bond:
            molList.append("%3d%3d%3d%3d  0  0  0\n" % (b.index[0]+1, 
                b.index[1]+1, b.order,b.stereo))

            # if necessary, write M  CHG records for charged atoms
        charge_atoms = []
        charge_values = []
        for a in model.atom:
            if a.formal_charge != 0:
                charge_atoms.append(a)
        if len(charge_atoms):
            c = 0
            for a in model.atom:
                a.index = c
                c = c + 1
            while len(charge_atoms) != 0:
                chg_set = charge_atoms[0:8]
                charge_atoms = charge_atoms[8:]
                tline = "M  CHG%3d" % (len(chg_set))
                for i in chg_set:
                    tline = tline + "%4d%4d" % (i.index+1,i.formal_charge)
                molList.append(tline + "\n")
        molList.append("M  END\n")
        return(molList)


