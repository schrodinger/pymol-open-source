from chempy.models import Indexed
from chempy import Storage,Atom,Bond

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
              gbasis='N31',ngauss=6,ndfunc=1):

         
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
      gmsList.append(" $DATA\n")
      gmsList.append(model.molecule.title+" 6-31G* optimization\n")
      gmsList.append("C1\n")

      # write atom records
      for a in model.atom:
         gmsList.append("%4s %5.2f %12.6f %12.6f %12.6f\n" %
                        (a.name,atNum[a.symbol],a.coord[0],
                         a.coord[1],a.coord[2]))

      gmsList.append(" $END\n")
      gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE. $END\n")
      gmsList.append(" $STATPT NPRT=-2 NPUN=-2 NSTEP=50 $END\n")
      gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
      gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
                     
      return(gmsList)







