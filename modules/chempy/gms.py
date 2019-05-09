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

#
# THIS MODULE IS OBSOLETE PLEASE DO NOT USE
#

import string

atNum = {
    'H'  : 1,
    'C'  : 6,
    'N'  : 7,
    'O'  : 8,
    'F'  : 9,
    'P'  : 15,
    'S'  : 16,
    'Cl' : 17,
    'Br' : 35,
    'I'  : 53,
    }

class GMS(Storage):

    def toList(self,model,runtyp='OPTIMIZE',exetyp='RUN',
                  gbasis='N31',ngauss=6,ndfunc=1,dirscf=1):


        gmsList = []

        # write header records

        nzvar = (model.nAtom*3)-6
        chg = 0
        for a in model.atom:
            chg = chg + a.formal_charge
        chg = int(chg)
        if chg==0:
            icharg=''
            diffsp=''
        else:
            icharg='ICHARG=%d' % chg
            if chg<0:
                diffsp='DIFFSP=.TRUE.'
            else:
                diffsp=''
        gmsList.append(
            " $CONTRL RUNTYP=%s COORD=UNIQUE EXETYP=%s NZVAR=%d %s $END\n" %
            (runtyp,exetyp,nzvar,icharg))
        if ndfunc>0:
            gmsList.append(" $BASIS GBASIS=%s NGAUSS=%d NDFUNC=%d %s $END\n" %
                                (gbasis,ngauss,ndfunc,diffsp))
        else:
            gmsList.append(" $BASIS GBASIS=%s NGAUSS=%d %s $END\n",
                                (gbasis,ngauss,diffsp))
        if dirscf:
            gmsList.append(" $SCF DIRSCF=.TRUE. $END\n")
        gmsList.append(" $DATA\n")
        gmsList.append(model.molecule.title+" 6-31G* optimization\n")
        gmsList.append("C1\n")

        # write atom records in an ordering compatible with internal
        # coordinate generation
        c = 1
        for z in model.get_internal_tuples():
            a = model.atom[z[0]]
            if not len(a.name):
                name = a.symbol + "%02d"%c
            else:
                name = a.name
            gmsList.append("%4s %5.2f %12.6f %12.6f %12.6f\n" %
                                (name,atNum[a.symbol],a.coord[0],
                                 a.coord[1],a.coord[2]))
            c = c + 1
        gmsList.append(" $END\n")
        gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE. $END\n")
        if runtyp=='OPTIMIZE':
            gmsList.append(" $STATPT NPRT=-2 NPUN=-2 NSTEP=50 $END\n")
        gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
        gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")

        return(gmsList)
