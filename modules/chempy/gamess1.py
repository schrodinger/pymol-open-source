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

# this is quick and dirty irst crack at a gamess interface, written by someone who
# knows very little about the program (WLD) - hence the version 1 identifier...

# by the way, most of the below is untested...

import os
import shutil
import glob
import re
import sys
import time

from chempy import feedback
from chempy.brick import Brick

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

def do(input,run_prefix=None,echo=None,
         punch=None,output=None,skip=None):
    if not run_prefix:
        run_prefix = 'gamess_run'
    if not skip:
        if feedback['gamess']:
            print(" "+str(__name__)+': creating temporary files "%s.*"' % (run_prefix))
            print(" "+str(__name__)+': launching gamess...')
        try:
            for a in glob.glob(run_prefix+".*"):
                os.unlink(a)
        except:
            pass
        f = open(run_prefix+".inp",'w')
        for a in input:
            f.write(a)
        f.close()
        if echo:
            os.system(rungms_path+' '+run_prefix+" 2>&1 | tee "+run_prefix+".out")
        else:
            os.system(rungms_path+' '+run_prefix+" > "+run_prefix+".out 2>&1")
# NFS workaround (flushes the directory cache so that glob will work)
        try: os.unlink(".sync")
        except: pass
        f = open(".sync",'w')
        f.close()
#
    if feedback['gamess']:
        print(" "+str(__name__)+': job complete. ')
    if punch:
        for src in glob.glob(run_prefix+".dat"):
            f = open(src)
            punch = f.readlines()
            f.close()
    if output:
        for src in glob.glob(run_prefix+".out"):
            f = open(src)
            output = f.readlines()
            f.close()
    return (output,punch)

if 'GAMESS' in os.environ:
    base = os.environ['GAMESS']
    bin_path = base + '/bin/'
    rungms_path = bin_path + 'rungms'
else:
    base = ''
    bin_path = ''
    params_path = ''

class State:

    def __init__(self):
        self.model = None
        self.data = None
        self.vec = None

    def load_model(self,model):
        self.model = model

    def get_zmat_ordering(self):
        lst = []
        for z in self.model.get_internal_tuples():
            lst.append(z[0])
        return lst

    def get_data_group(self,basis = None,zmat = 1):
        model = self.model

        gmsList = []

        # write header records
        gmsList.append(" $DATA\n")
        gmsList.append(model.molecule.title+" from "+str(__name__)+"\n")
        gmsList.append("C1\n")

        # write atom records in an ordering compatible with internal
        # coordinate generation
        c = 1
        for z in self.get_zmat_ordering():
            a = model.atom[z]
            if not len(a.name):
                name = a.symbol + "%02d"%c
            else:
                name = a.name
            gmsList.append("%10s %5.1f %18.10f %18.10f %18.10f\n" %
                                (name,atNum[a.symbol],a.coord[0],
                                 a.coord[1],a.coord[2]))
            c = c + 1
        gmsList.append(" $END\n")
        return gmsList

    def get_ordered_data_group(self):
        gmsList = self.data[0:3]
        flag = 1
        c = 3
        for a in self.data[3:]:
            if flag:
                flag = 0
                gmsList.append(a)
            if not a.strip():
                flag = 1
            c = c + 1
        return gmsList

    def get_contrl_group(self,
                                scftyp='RHF',
                                runtyp='ENERGY',
                                exetyp='RUN',
                                coord='UNIQUE',
                                nzvar = -1):
        gmsList = []
        model = self.model
        if nzvar:
            if nzvar<0:
                nzvar = (self.model.nAtom*3)-6
        gmsList.append(" $CONTRL SCFTYP=%s RUNTYP=%s EXETYP=%s\n"
                            % (scftyp,runtyp,exetyp) )
        if coord:
            gmsList.append("COORD=%s\n"%coord)
        if nzvar:
            gmsList.append("NZVAR=%d\n"%nzvar)
        chg = 0
        for a in model.atom:
            chg = chg + a.formal_charge
        chg = int(chg)
        if chg==0:
            icharg=None
        else:
            icharg='ICHARG=%d' % chg
        if icharg:
            gmsList.append("%s\n" % (icharg))
        gmsList.append(" $END\n")
        return gmsList

    def read_output_list(self,list):
        ll = len(list)
        c = 0
        crd_list = []
        chg_list = []
        nrg_list = []
        for a in list:
            if a[0:36] == ' COORDINATES OF ALL ATOMS ARE (ANGS)':
                crd_list.append(c+3)
            if a[0:13] == ' NET CHARGES:':
                chg_list.append(c+4)
            if a[0:37] == '                       TOTAL ENERGY =':
                nrg_list.append(c)
            c = c + 1
        atom = self.model.atom
        idx = {}
        c = 0
        for a in atom:
            idx[a.name.upper()]=c  # games converts to uppercase
            c = c + 1
        if len(crd_list):
            a = crd_list.pop()
            cc = 0
            while a<ll:
                l = list[a]
                name = l[1:11].strip()
                if name=='':
                    break
                atom[idx[name]].coord = [float(l[16:31]),
                                                 float(l[31:46]),
                                                 float(l[46:61])]
                cc = cc + 1
                a = a + 1
            if cc and feedback['gamess']:
                print(" "+str(__name__)+': coordinates modified for %d atoms.' % (cc))
        if len(chg_list):
            a = chg_list.pop()
            cc = 0
            while a<ll:
                l = list[a]
                name = l[1:11].strip()
                if name[0]=='-':
                    break
                atom[idx[name]].partial_charge = float(l[19:27])
                a = a + 1
                cc = cc + 1
            if cc and feedback['gamess']:
                print(" "+str(__name__)+': charges modified for %d atoms.' % (cc))
        if len(nrg_list):
            a = nrg_list.pop()
            l = list[a]
            # get energy, and convert to kcal/mole
            self.model.molecule.energy = float(l[38:58].strip())*627.5095
            if feedback['gamess']:
                print(" "+str(__name__)+': energy updated %12.6f.' % self.model.molecule.energy)

    def read_punch_list(self,list):
        ll = len(list)
        c = 0
        data_list = []
        vec_list = []
        for a in list:
            if a[0:6] == ' $DATA':
                data_list.append(c)
            elif a[0:5] == ' $VEC':
                vec_list.append(c)
            c = c + 1
        if len(data_list):
            a = data_list.pop()
            self.data = []
            data = self.data
            while a<ll:
                la = list[a]
                data.append(la)
                if la[0:5] == ' $END':
                    break
                a = a + 1
            if feedback['gamess']:
                print(" "+str(__name__)+': read $DATA group.')
        if len(vec_list):
            a = vec_list.pop()
            self.vec = []
            vec = self.data
            while a<ll:
                la = list[a]
                vec.append(la)
                if la[0:5] == ' $END':
                    break
                a = a + 1
            if feedback['gamess']:
                print(" "+str(__name__)+': read new $VEC group.')

    def update_data_coords(self): # update coordinates of ordered atoms in $DATA
        idx = {}
        c = 0
        for a in self.model.atom:
            idx[a.name.upper()]=c
            c = c + 1
        if self.data:
            flag = 1
            c = 3
            for a in self.data[3:]:
                if flag:
                    flag = 0
                    kee = a[0:3]
                    if kee not in idx:
                        break
                    i = idx[kee]
                    at = self.model.atom[i]
                    self.data[c]="%-10s%5.1f%18.10f%18.10f%18.10f\n" % (
                        at.name,atNum[at.symbol],at.coord[0],
                        at.coord[1],at.coord[2])
                if not a.strip():
                    flag = 1
                c = c + 1

    def read_density_list(self,list,brick,z_step):
        ll = len(list)
        c = 0
        den_list = []
        for a in list:
            if a[0:37] == ' ELECTRON DENSITY, IPOINT,X,Y,Z,EDENS':
                den_list.append(c+1)
            c = c + 1
        if len(den_list):
            lst = 0
            a = den_list.pop()
            for x in range(brick.dim[0]):
                for y in range(brick.dim[1]):
                    brick.lvl[x][y][z_step] = float(list[a][36:51])
                    a = a + 1
            if feedback['gamess']:
                print(" "+str(__name__)+': read density slice %d of %d.' %(
                    z_step+1,brick.dim[2]))

    def read_potential_list(self,list,brick,z_step):
        ll = len(list)
        c = 0
        pot_list = []
        for a in list:
            if a[0:51] == 'THE ROWS OF THE ELECTROSTATIC POTENTIAL GRID (A.U.)':
                pot_list.append(c+1)
            c = c + 1
        if len(pot_list):
            lst = 0
            a = pot_list.pop()
            for x in range(brick.dim[0]):
                aa = a
                mat = list[aa][0:3]
                col = []
                while 1:
                    if list[aa][0:3]==mat:
                        col.append(list[aa][5:])
                        aa = aa + 1
                    else:
                        break
                a = aa
                vst = ' '.join(col).split()
                for y in range(brick.dim[1]):
                    brick.lvl[x][y][z_step] = float(vst[y])
            if feedback['gamess']:
                print(" "+str(__name__)+': read potential slice %d of %d.' %(
                    z_step+1,brick.dim[2]))

    def get_basis_group(self,gbasis='N31',ngauss=6,ndfunc=1):
        gmsList = []
        gmsList.append(" $BASIS GBASIS=%s NGAUSS=%d NDFUNC=%d\n" %
                            (gbasis,ngauss,ndfunc))
        model = self.model
        chg = 0
        for a in model.atom:
            chg = chg + a.formal_charge
        chg = int(chg)
        if chg<0:
            diffsp=' DIFFSP=.TRUE.'
        else:
            diffsp=None
        if diffsp:
            gmsList.append("%s\n" % (diffsp))
        gmsList.append(" $END\n")
        return gmsList

    def get_zmat_ext_freeze_torsion(self,flag=3):
        # requires PYMOL to read dihedrals from structure
        # requires list of dihedrals from tinker.amber
        #
        from pymol import cmd
        from .tinker.amber import Topology

        cmd.load_model(self.model,'_gamess1')
        model = self.model

        # get mapping of model ordering to zmat ordering
        m2z = {}
        z2m = {}
        c = 1 # GAMESS is one-based
        for a in self.get_zmat_ordering():
            m2z[a] = c
            z2m[c] = a
            c = c + 1

        # get all torsions in the molecule

        topo = Topology(self.model)

        # find those where flag is set in all atoms

        mask = 2 ** flag

        frozen_list = []

        for a in list(topo.torsion.keys()):
            if (model.atom[a[0]].flags&
                 model.atom[a[1]].flags&
                 model.atom[a[2]].flags&
                 model.atom[a[3]].flags)&mask:
                frozen_list.append(a)

        print(" freeze-torsion: %d torsions will be frozen."%len(frozen_list))

        irzmat = []
        ifzmat = []
        fvalue = []
        if len(frozen_list):

            for frozen in frozen_list:
                # find additional torsions which need to be removed

                remove = []

                for a in list(topo.torsion.keys()):
                    if (((a[1]==frozen[1])and(a[2]==frozen[2])) or
                         ((a[2]==frozen[1])and(a[1]==frozen[2]))):
                        if a!=frozen:
                            remove.append(a)

                # convert to internal coordinate ordering

                frozen_z = (m2z[frozen[0]],m2z[frozen[1]],
                                m2z[frozen[2]],m2z[frozen[3]])

                remove_z = []
                for a in remove:
                    remove_z.append(m2z[a[0]],m2z[a[1]],m2z[a[2]],m2z[a[3]])

                # now reorder atoms in torsions to reflect z_matrix ordering
                # (not sure this is necessary)

                if frozen_z[0]>frozen_z[3]:
                    frozen_z = (frozen_z[3],frozen_z[2],frozen_z[1],frozen_z[0])
                tmp_z = []
                for a in remove_z:
                    if a[0]>a[3]:
                        tmp_z.append((a[3],a[2],a[1],a[0]))
                    else:
                        tmp_z.append(a)
                remove_z = tmp_z

                # get value of the fixed torsion

                fixed = (z2m[frozen_z[0]],z2m[frozen_z[1]],
                            z2m[frozen_z[2]],z2m[frozen_z[3]])

                dihe = cmd.get_dihedral("(_gamess1 and id %d)"%fixed[0],
                                      "(_gamess1 and id %d)"%fixed[1],
                                      "(_gamess1 and id %d)"%fixed[2],
                                      "(_gamess1 and id %d)"%fixed[3])

                # write out report for user edification

                print(" freeze-torsion: fixing freeze-torsion:")
                print(" freeze-torsion: %d-%d-%d-%d (pymol), %d-%d-%d-%d (gamess)"%(
                    fixed[0],fixed[1],fixed[2],fixed[3],
                    frozen_z[0],frozen_z[1],frozen_z[2],frozen_z[3]))
                print(" freeze-torsion: at %5.3f"%dihe)
                print(" freeze-torsion: removing redundant torsions:")
                for a in remove_z[1:]:
                    print(" freeze-torsion: %d-%d-%d-%d (pymol), %d-%d-%d-%d (gamess)"%(
                        z2m[a[0]],z2m[a[1]],z2m[a[2]],z2m[a[3]],
                        a[0],a[1],a[2],a[3]))

                # add parameters for this torsion into the list

                ifzmat.append((3,frozen_z[0],frozen_z[1],frozen_z[2],frozen_z[3]))
                fvalue.append(dihe)

                if len(remove_z):
                    for a in remove_z[1:]:
                        irzmat.append((3,a[0],a[1],a[2],a[3]))

        # generate restrained dihedral information

        zmat_ext = []
        if len(ifzmat):
            zmat_ext.append(" IFZMAT(1)=\n")
            comma = ""
            for a in ifzmat:
                zmat_ext.append(comma+"%d,%d,%d,%d,%d\n"%a)
                comma = ","
        if len(fvalue):
            zmat_ext.append(" FVALUE(1)=\n")
            comma = ""
            for a in fvalue:
                zmat_ext.append(comma+"%1.7f\n"%a)
                comma = ","
        if len(irzmat):
            zmat_ext.append(" IRZMAT(1)=\n")
            comma = ""
            for a in irzmat:
                zmat_ext.append(comma+"%d,%d,%d,%d,%d\n"%a)
                comma = ","

        cmd.delete("_gamess1") # important
        if len(zmat_ext):
            return zmat_ext
        else:
            return None

    def get_zmat_group(self,auto=1,dlc=1,zmat_extend=None):
        gmsList = []
        if auto and dlc:
            if zmat_extend is None:
                gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE. $END\n")
            else:
                gmsList.append(" $ZMAT DLC=.TRUE. AUTO=.TRUE.\n")
                gmsList.extend(zmat_extend)
                gmsList.append(" $END\n")
        else:
            raise RuntimeError
        return gmsList

    def get_eldens_group(self,morb=0):
        gmsList = []
        gmsList.append(" $ELDENS IEDEN=1 MORB=%i \n" % morb)
        gmsList.append("WHERE=GRID OUTPUT=PUNCH\n $END\n")
        return gmsList

    def get_elpot_group(self,morb=0):
        gmsList = []
        gmsList.append(" $ELPOT IEPOT=1 \n")
        gmsList.append("WHERE=GRID OUTPUT=PUNCH\n $END\n")
        return gmsList

    def get_guess_group(self,guess='HUCKEL'):
        return [" $GUESS GUESS=%s $END\n"%guess]

    def get_grid_group(self,brick,z_step):
        origin = (
            brick.origin[0],
            brick.origin[1],
            brick.origin[2]+brick.grid[2]*z_step)
        x_coord = (
            brick.origin[0]+brick.range[0]+brick.grid[0]/100.0,
            brick.origin[1],
            brick.origin[2]+brick.grid[2]*z_step)
        y_coord = (
            brick.origin[0],
            brick.origin[1]+brick.range[1]+brick.grid[1]/100.0,
            brick.origin[2]+brick.grid[2]*z_step)
        gmsList = [
            " $GRID ORIGIN(1)=%12.5f,%12.5f,%12.5f\n" % origin,
            "XVEC(1) = %12.5f,%12.5f,%12.5f\n" % x_coord,
            "YVEC(1) = %12.5f,%12.5f,%12.5f\n" % y_coord,
            "SIZE = %12.5f\n" % brick.grid[0],
            " $END\n"
            ]
        return gmsList

    def get_scf(self,dirscf=1):
        gmsList = []
        if dirscf:
            gmsList.append(" $SCF DIRSCF=.TRUE. $END\n")
        return gmsList

    def get_optimize_job(self,dirscf=1,zmat_extend=None):
        gmsList = []
        gmsList.extend(self.get_contrl_group(runtyp='OPTIMIZE'))
        gmsList.extend(self.get_basis_group())
        gmsList.extend(self.get_scf(dirscf=dirscf))
        gmsList.extend(self.get_data_group())
        gmsList.extend(self.get_zmat_group(zmat_extend=zmat_extend))
        gmsList.append(" $STATPT NSTEP=50 $END\n")
        return gmsList

    def get_optimize_charge_job(self):
        gmsList = self.get_optimize_job()
        gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
        gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
        return gmsList

    def get_energy_charge_job(self):
        gmsList = self.get_energy_job()
        gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
        gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
        return gmsList

    def get_prop_job(self):
        gmsList = []
        gmsList.extend(self.get_contrl_group(runtyp = 'PROP',
                                                         nzvar=0))
        gmsList.extend(self.get_guess_group(guess='MOREAD'))
        self.update_data_coords()
        gmsList.extend(self.data)
        gmsList.extend(self.vec)
        return gmsList

    def get_energy_job(self):
        gmsList=[]
        gmsList.extend(self.get_contrl_group(
            runtyp = 'ENERGY'
            ))
        gmsList.extend(self.get_basis_group())
        gmsList.extend(self.get_scf())
        gmsList.extend(self.get_data_group())
        gmsList.extend(self.get_zmat_group())
        return gmsList

    def get_density_job(self,brick,z_step,morb=0):
        gmsList = self.get_prop_job()
        gmsList.extend(self.get_eldens_group(morb=morb))
        gmsList.extend(self.get_grid_group(brick,z_step))
        return gmsList

    def get_potential_job(self,brick,z_step):
        gmsList = self.get_prop_job()
        gmsList.extend(self.get_elpot_group())
        gmsList.extend(self.get_grid_group(brick,z_step))
        return gmsList

    def get_prop_charge_job(self):
        gmsList = self.get_prop_job()
        gmsList.append(" $ELPOT IEPOT=1 WHERE=PDC $END\n")
        gmsList.append(" $PDC PTSEL=GEODESIC CONSTR=CHARGE $END\n")
        return gmsList

    def get_density(self,brick,morb=0,run_prefix=None):
        for a in range(brick.dim[2]):
            gmsList = self.get_density_job(brick,a,morb=morb)
            result = do(gmsList,punch=1,
                            run_prefix=run_prefix)
            self.read_density_list(result[1],brick,a)

    def get_potential(self,brick,run_prefix=None):
        for a in range(brick.dim[2]):
            gmsList = self.get_potential_job(brick,a)
            result = do(gmsList,punch=1,
                            run_prefix=run_prefix)
            self.read_potential_list(result[1],brick,a)

    def get_charges(self,run_prefix=None):
        gmsList = self.get_energy_job()
        result = do(gmsList,output=1,punch=1,
                        run_prefix=run_prefix)
        self.read_output_list(result[0])
        self.read_punch_list(result[1])

    def get_energy(self,run_prefix=None):
        gmsList = self.get_energy_job()
        result = do(gmsList,output=1,punch=1,
                        run_prefix=run_prefix)
        self.read_output_list(result[0])
        self.read_punch_list(result[1])

    def get_optimized_energy(self,run_prefix=None,zmat_extend=None):
        gmsList = self.get_optimize_job(zmat_extend=zmat_extend)
        result = do(gmsList,output=1,punch=1,
                        run_prefix=run_prefix)
        self.read_output_list(result[0])
        self.read_punch_list(result[1])

    def get_optimized_charges(self,run_prefix=None,skip=None):
        gmsList = self.get_optimize_charge_job()
        result = do(gmsList,output=1,punch=1,
                        run_prefix=run_prefix,skip=skip)
        self.read_output_list(result[0])
        self.read_punch_list(result[1])

    def get_prop_charges(self,run_prefix=None):
        gmsList = self.get_prop_charge_job()
        result = do(gmsList,output=1,punch=1,
                        run_prefix=run_prefix)
        self.read_output_list(result[0])
        self.read_punch_list(result[1])


if 'GAMESS' in os.environ:
    base = os.environ['GAMESS']
    rungms_path = base + '/bin/rungms'
else:
    base = ''
    rungms_path = ''
