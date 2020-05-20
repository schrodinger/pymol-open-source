from pymol.wizard import Wizard
from pymol import cmd,editor
from chempy import io
from copy import deepcopy

import pymol
import os
import traceback

src_sele = "_mutate_sel"
bump_name = "_bump_check"

obj_name = "mutation"
frag_name = "_tmp_mut"
mut_sele = "_tmp_mut_sele"
tmp_obj2 = "_tmp_obj2"
tmp_sele1 = "_tmp_sele1"
tmp_sele2 = "_tmp_sele2"
tmp_hbonds = "_tmp_hbonds"

default_mode = "current"
default_rep = "lines"
default_hyd = 'auto'
default_dep = 'dep'
default_n_cap = 'none'
default_c_cap = 'none'

_rot_type_xref = {
    'GLUH' : 'GLU',
    'ASPH' : 'ASP',
    'ARGN' : 'ARG',
    'LYSN' : 'LYS',
    'HIP' : 'HIS',
    'HID' : 'HIS',
    'HIP' : 'HIS'
    }

class Mutagenesis(Wizard):

    count = 0
    cutoff = 3.5

    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)
        cmd=self.cmd

        if self.cmd.get_movie_length() > 0:
            raise pymol.wizarding.WizardError('Mutagenesis Wizard cannot be used with Movie')

        cmd.unpick()

        self.stored = pymol.Scratch_Storage()
        self.space = {'stored': self.stored}

        self.bump_scores = []
        self.dep = default_dep

        self.ind_library = io.pkl.fromFile(os.environ['PYMOL_DATA']+
                                           "/chempy/sidechains/sc_bb_ind.pkl")
        self.load_library()
        self.status = 0 # 0 no selection, 1 mutagenizing
        self.bump_check = 1
        self.auto_center = 1
        self.error = None
        self.object_name = None
        self.modes = [
            'current'
            ]
        self.mode = default_mode
        self.rep = default_rep
        self.hyd = default_hyd
        self.n_cap = default_n_cap
        self.c_cap = default_c_cap
        residues = list(self.ind_library.keys())
        # could extent with additional fragments manually as below
        residues.extend(['GLY','ALA'])
        residues.extend(['HID','HIE','HIP'])
        residues.extend(['ARGN','LYSN','ASPH','GLUH'])
        residues.sort()
        res_copy = deepcopy(residues)
        for a in res_copy:
            residues.append('NT_'+a)
            residues.append('CT_'+a)
        self.modes.extend(residues)
        self.mode_label={}
        for a in self.modes:
            self.mode_label[a] = ""+a
        self.mode_label['current']="No Mutant"

        self.selection_mode = cmd.get_setting_int("mouse_selection_mode")
        cmd.set("mouse_selection_mode",1)

        smm = []
        smm.append([ 2, 'Mutant', '' ])
        smm.append([ 1, 'No change', 'cmd.get_wizard().set_mode("current")' ])
#        smm.append([ 1, 'N-Term', [] ])
#        smm.append([ 1, 'C-Term', [] ])
        smm.append([ 0, '', '' ])
        for a in self.modes:
            if a == 'current':
                pass
            elif a[0:3]=='NT_':
                pass
#                smm[2][2].append([ 1, self.mode_label[a[3:]], 'cmd.get_wizard().set_mode("'+a+'")'])
            elif a[0:3]=='CT_':
                pass
#                smm[3][2].append([ 1, self.mode_label[a[3:]], 'cmd.get_wizard().set_mode("'+a+'")'])
            else:
                smm.append([ 1, self.mode_label[a], 'cmd.get_wizard().set_mode("'+a+'")'])

        # group arg, lys, his, glu, asp

        for lst in [ smm ]: # [ smm, smm[2][2], smm[3][2] ]:
            for a in 'ARG','LYS','HID','GLU','ASP':
                ix = 0
                start = 0
                stop = 0
                for b in lst:
                    if start==0:
                        if b[1][0:]==a:
                            start = ix
                            stop = ix + 1
                    elif b[1][0:3]==a[0:3] or ( b[1][0:2]==a[0:2] and a[0:2]=='HI' ):
                        stop = ix + 1
                    ix = ix + 1
                if start!=0 and stop!=0:
                    slice = lst[start:stop]
                    if a != 'HID':
                        slice2 = [slice[0] ] + [ [0,'',''] ] + slice[1:]
                        lst[start:stop] = [ [1, self.mode_label[a] + "... " , slice2 ] ]
                    else:
                        slice2 = [ slice[3] ] + [ [0,'',''] ] + slice[0:3]
                        lst[start:stop] = [ [1, self.mode_label['HIS']+ "... ", slice2 ] ]

        self.menu['mode']=smm


        self.reps = [
            'lines',
            'sticks',
            'spheres',
            'dots'
            ]

        self.rep_name = {
            'lines' : "Show Lines",
            'sticks' : "Show Sticks",
            'spheres' : "Show Spheres",
            'dots' : "Show Dots",
            }

        self.dep_name = {
            'dep' : "Backbone Depen. Rotamers",
            'ind' : "Backbone Indep. Rotamers"
            }

        self.hyd_name = {
            'auto' : "Hydrogens: Current",
            'keep' : "Hydrogens: Add & Retain",
#            'polar' : "Polar Hydrogens",
            'none'  : "Hydrogens: Remove",
            }
        self.hyds = [ 'auto', 'keep', 'none' ]

        self.n_cap_name = {
            'none' : 'Open',
            'posi' : 'NH3+',
            'acet' : 'Acetyl',
            }
        self.n_caps = [ 'none', 'posi', 'acet' ]

        self.c_cap_name = {
           'none' : 'Open',
           'nega' : 'COO-',
           'amin' : 'Amine',
           'nmet' : 'N-methyl',
            }
        self.c_caps = [ 'none', 'nega', 'amin', 'nmet' ]

        smm = []
        smm.append([ 2, 'N-Cap', '' ])
        for a in self.n_caps:
            smm.append([ 1, self.n_cap_name[a], 'cmd.get_wizard().set_n_cap("'+a+'")'])
        self.menu['n_cap']=smm

        smm = []
        smm.append([ 2, 'C-Cap', '' ])
        for a in self.c_caps:
            smm.append([ 1, self.c_cap_name[a], 'cmd.get_wizard().set_c_cap("'+a+'")'])
        self.menu['c_cap']=smm

        smm = []
        smm.append([ 2, 'Hydrogens', '' ])
        for a in self.hyds:
            smm.append([ 1, self.hyd_name[a], 'cmd.get_wizard().set_hyd("'+a+'")'])
        self.menu['hyd']=smm

        smm = []
        smm.append([ 2, 'Representation', '' ])
        for a in self.reps:
            smm.append([ 1, self.rep_name[a], 'cmd.get_wizard().set_rep("'+a+'")'])
        self.menu['rep']=smm

        self.deps = [ 'dep', 'ind' ]
        smm = []
        smm.append([ 2, 'Rotamers', '' ])
        for a in self.deps:
            smm.append([ 1, self.dep_name[a], 'cmd.get_wizard().set_dep("'+a+'")'])
        self.menu['dep']=smm

        if 'pk1' in cmd.get_names('selections'):
            cmd.select(src_sele,"(byres pk1)")
            cmd.unpick()
            cmd.enable(src_sele)
            self.status = 1
            self.error = None
            self.do_library()
            cmd.refresh_wizard()

    def load_library(self):
        if self.dep == 'dep':
            if not hasattr(self,'dep_library'):
                self.dep_library = io.pkl.fromFile(os.environ['PYMOL_DATA']+
                                           "/chempy/sidechains/sc_bb_dep.pkl")

    def set_mode(self,mode):
        cmd=self.cmd
        if mode in self.modes:
            self.mode = mode
        if self.status==1:
            self.do_library()
        cmd.refresh_wizard()

    def set_rep(self,rep):
        cmd=self.cmd
        if rep in self.reps:
            self.rep=rep
        cmd.hide("("+obj_name+")")
        cmd.show('lines',obj_name) # always show lines
        cmd.show(self.rep,obj_name)
        cmd.refresh_wizard()

    def set_c_cap(self,c_cap):
        cmd=self.cmd
        if c_cap in self.c_caps:
            self.c_cap=c_cap
        if self.status==1:
            self.do_library()
        cmd.refresh_wizard()

    def set_n_cap(self,n_cap):
        cmd=self.cmd
        if n_cap in self.n_caps:
            self.n_cap=n_cap
        if self.status==1:
            self.do_library()
        cmd.refresh_wizard()

    def set_hyd(self,hyd):
        cmd=self.cmd
        if hyd in self.hyds:
            self.hyd=hyd
        if self.status==1:
            self.do_library()
        cmd.refresh_wizard()

    def set_dep(self,value):
        cmd=self.cmd
        if value!=self.dep:
            self.dep = value
            self.load_library()
            if src_sele in cmd.get_names("all"):
                self.do_library()
            cmd.refresh_wizard()

    def get_panel(self):
        cmd=self.cmd
        if int(cmd.get("mouse_selection_mode")!=1):
            cmd.set("mouse_selection_mode",1)
        if self.mode == 'current':
            label = 'No Mutation'
        else:
            label = 'Mutate to '+self.mode_label[self.mode]
        return [
            [ 1, 'Mutagenesis',''],
            [ 3, label,'mode'],
            [ 3, 'N-Cap: '+self.n_cap_name[self.n_cap],'n_cap'],
            [ 3, 'C-Cap: '+self.c_cap_name[self.c_cap],'c_cap'],
            [ 3, self.hyd_name[self.hyd],'hyd'],
            [ 3, self.rep_name[self.rep],'rep'],
            [ 3, self.dep_name[self.dep],'dep'],
            [ 2, 'Apply' , 'cmd.get_wizard().apply()'],
            [ 2, 'Clear' , 'cmd.get_wizard().clear()'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def get_event_mask(self):
        return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_state

    def cleanup(self):
        cmd=self.cmd
        global default_mode,default_rep,default_dep,default_hyd
        global default_n_cap, default_c_cap
        default_mode = self.mode
        default_rep = self.rep
        default_dep = self.dep
        default_hyd = self.hyd
        default_n_cap = self.n_cap
        default_c_cap = self.c_cap
        cmd.set("mouse_selection_mode",self.selection_mode) # restore selection mode
        self.clear()

    def clear(self):
        cmd=self.cmd
        self.status=0
        self.bump_scores = []
        cmd.delete(tmp_hbonds)
        cmd.delete(tmp_obj2)
        cmd.delete(mut_sele)
        cmd.delete(src_sele)
        cmd.delete(obj_name)
        cmd.delete(bump_name)
        cmd.delete("_seeker_hilight")
        cmd.refresh_wizard()

    def apply(self):
        cmd=self.cmd
        if self.status==1:
            # find the name of the object which contains the selection
            src_frame = cmd.get_state()
            try:
                new_name = cmd.get_object_list(src_sele)[0]
            except IndexError:
                print(" Mutagenesis: object not found.")
                return

            if True:
                auto_zoom = cmd.get_setting_text('auto_zoom')
                cmd.set('auto_zoom',"0",quiet=1)
                if self.lib_mode!="current":
                    # create copy with mutant in correct frame
                    state = cmd.get_object_state(new_name)
                    cmd.create(tmp_obj2, obj_name, src_frame, state)
                    cmd.set_title(tmp_obj2, state, '')
                    cmd.color(self.stored.identifiers[4], "?%s & elem C" % tmp_obj2)
                    cmd.alter(tmp_obj2, 'ID = -1')

                    # select backbone connection atoms
                    cmd.select(tmp_sele1, 'neighbor ?%s' % (src_sele), 0)

                    # remove residue and neighboring c-cap/n-cap (if any)
                    cmd.remove("?%s | byres (?%s & "
                            "(name N & resn NME+NHH | name C & resn ACE))" % (src_sele, tmp_sele1))

                    # create the merged molecule
                    cmd.create(new_name, "?%s | ?%s" % (new_name, tmp_obj2), state, state)

                    # now connect them
                    cmd.select(tmp_sele2, '/%s/%s/%s/%s' % ((new_name,) + self.stored.identifiers[:3]))
                    cmd.bond('?%s & name C' % (tmp_sele1), '?%s & name N' % (tmp_sele2), quiet=1)
                    cmd.bond('?%s & name N' % (tmp_sele1), '?%s & name C' % (tmp_sele2), quiet=1)
                    cmd.set_geometry('(?%s | ?%s) & name C+N' % (tmp_sele1, tmp_sele2), 3, 3) # make amide planer

                    # fix N-H hydrogen position (if any exists)
                    cmd.h_fix('?%s & name N' % (tmp_sele2))

                    # delete temporary objects/selections
                    cmd.delete(tmp_sele1)
                    cmd.delete(tmp_sele2)
                    cmd.delete(tmp_obj2)
                    self.clear()
                    # and return to frame 1
                    cmd.frame(1)
                    cmd.refresh_wizard()
                else:
                    # create copy with conformation in correct state
                    cmd.create(tmp_obj2,obj_name,src_frame,1)

                    # remove existing c-cap in copy (if any)
                    cmd.remove("byres (name N and (%s in (neighbor %s)) and resn NME+NHH)"%
                                (new_name,src_sele))
                    cmd.remove("(%s) and name OXT"%src_sele)

                    # remove existing n-cap in copy (if any)
                    cmd.remove("byres (name C and (%s in (neighbor %s)) and resn ACE)"%
                                (new_name,src_sele))

                    # save existing conformation on undo stack
#               cmd.edit("((%s in %s) and name ca)"%(new_name,src_sele))
                    cmd.push_undo("("+src_sele+")")
                    # modify the conformation
                    cmd.update(new_name,tmp_obj2)
#               cmd.unpick()
                    cmd.delete(tmp_obj2)
                    self.clear()
                    # and return to frame 1
                    cmd.frame(1)
                    cmd.refresh_wizard()
                cmd.set('auto_zoom',auto_zoom,quiet=1)

    def get_prompt(self):
        self.prompt = None
        if self.status==0:
            self.prompt = [ 'Pick a residue...']
        elif self.status==1:
            self.prompt = [ 'Select a rotamer for %s or pick a new residue...'%self.res_text ]
        return self.prompt


    def do_library(self):
        cmd=self.cmd
        pymol=cmd._pymol
        if not ((cmd.count_atoms("(%s) and name N"%src_sele)==1) and
                (cmd.count_atoms("(%s) and name C"%src_sele)==1) and
                (cmd.count_atoms("(%s) and name O"%src_sele)==1)):
            self.clear()
            return 1
        cmd.feedback("push")
        cmd.feedback("disable","selector","everythin")
        cmd.feedback("disable","editor","actions")
        self.prompt = [ 'Loading rotamers...']
        self.bump_scores = []
        state_best = 0

        pymol.stored.name = 'residue'
        cmd.iterate("first (%s)"%src_sele,'stored.name=model+"/"+segi+"/"+chain+"/"+resn+"`"+resi')
        self.res_text = pymol.stored.name
        cmd.select("_seeker_hilight",src_sele)

        auto_zoom = cmd.get_setting_text('auto_zoom')
        cmd.set('auto_zoom',"0",quiet=1)
        cmd.frame(0)
        cmd.delete(frag_name)
        if self.auto_center:
            cmd.center(src_sele,animate=-1)

        self.lib_mode = self.mode
        if self.lib_mode=="current":
            pymol.stored.resn=""
            cmd.iterate("(%s & name CA)"%src_sele,"stored.resn=resn")
            rot_type = _rot_type_xref.get(pymol.stored.resn,pymol.stored.resn)
            if (self.c_cap!='none') or (self.n_cap!='none') or (self.hyd != 'auto'):
                self.lib_mode = rot_type # force fragment-based load
            else:
                cmd.create(frag_name,src_sele,1,1)
                if self.c_cap=='open':
                    cmd.remove("%s and name OXT"%frag_name)

        if self.lib_mode!='current':
            rot_type = self.lib_mode
            frag_type = self.lib_mode
            if (self.n_cap == 'posi') and (frag_type[0:3]!='NT_'):
                if not ( cmd.count_atoms(
                    "elem C & !(%s) & (bto. (name N & (%s))) &! resn ACE"%
                                     (src_sele,src_sele))):
                    # use N-terminal fragment
                    frag_type ="NT_"+frag_type
            if (self.c_cap == 'nega') and (frag_type[0:3]!='CT_'):
                if not ( cmd.count_atoms("elem N & !(%s) & (bto. (name C & (%s))) & !resn NME+NHH"%
                                     (src_sele,src_sele))):
                    # use C-terminal fragment
                    frag_type ="CT_"+frag_type
            if rot_type[0:3] in [ 'NT_', 'CT_' ]:
                rot_type = rot_type[3:]
            rot_type = _rot_type_xref.get(rot_type, rot_type)
            cmd.fragment(frag_type.lower(), frag_name, origin=0)
            # trim off hydrogens
            if (self.hyd == 'none'):
                cmd.remove("("+frag_name+" and hydro)")
            elif (self.hyd == 'auto'):
                if cmd.count_atoms("("+src_sele+") and hydro")==0:
                    cmd.remove("("+frag_name+" and hydro)")
            # copy identifying information
            cmd.alter("?%s & name CA" % src_sele, "stored.identifiers = (segi, chain, resi, ss, color)", space=self.space)
            cmd.alter("?%s" % frag_name, "(segi, chain, resi, ss) = stored.identifiers[:4]", space=self.space)
            # move the fragment
            if ((cmd.count_atoms("(%s & name CB)"%frag_name)==1) and
                 (cmd.count_atoms("(%s & name CB)"%src_sele)==1)):
                cmd.pair_fit("(%s & name CA)"%frag_name,
                             "(%s & name CA)"%src_sele,
                             "(%s & name CB)"%frag_name,
                             "(%s & name CB)"%src_sele,
                             "(%s & name C)"%frag_name,
                             "(%s & name C)"%src_sele,
                             "(%s & name N)"%frag_name,
                             "(%s & name N)"%src_sele)
            else:
                cmd.pair_fit("(%s & name CA)"%frag_name,
                             "(%s & name CA)"%src_sele,
                             "(%s & name C)"%frag_name,
                             "(%s & name C)"%src_sele,
                             "(%s & name N)"%frag_name,
                             "(%s & name N)"%src_sele)

            # fix the carbonyl position...
            cmd.iterate_state(1,"(%s & name O)"%src_sele,"stored.list=[x,y,z]")
            cmd.alter_state(1,"(%s & name O)"%frag_name,"(x,y,z)=stored.list")
            if cmd.count_atoms("(%s & name OXT)"%src_sele):
                cmd.iterate_state(1,"(%s & name OXT)"%src_sele,"stored.list=[x,y,z]")
                cmd.alter_state(1,"(%s & name OXT)"%frag_name,"(x,y,z)=stored.list")
            elif cmd.count_atoms("(%s & name OXT)"%frag_name): # place OXT if no template exists
                angle = cmd.get_dihedral("(%s & name N)"%frag_name,
                                         "(%s & name CA)"%frag_name,
                                         "(%s & name C)"%frag_name,
                                         "(%s & name O)"%frag_name)
                cmd.protect("(%s & name O)"%frag_name)
                cmd.set_dihedral("(%s & name N)"%frag_name,
                                 "(%s & name CA)"%frag_name,
                                 "(%s & name C)"%frag_name,
                                 "(%s & name OXT)"%frag_name,180.0+angle)
                cmd.deprotect(frag_name)


            # fix the hydrogen position (if any)
            if cmd.count_atoms("(hydro and bound_to (name N & (%s)))"%frag_name)==1:
                if cmd.count_atoms("(hydro and bound_to (name N & (%s)))"%src_sele)==1:
                    cmd.iterate_state(1,"(hydro and bound_to (name N & (%s)))"%src_sele,
                                      "stored.list=[x,y,z]")
                    cmd.alter_state(1,"(hydro and bound_to (name N & (%s)))"%frag_name,
                                    "(x,y,z)=stored.list")
                elif cmd.select(tmp_sele1,"(name C & bound_to (%s and elem N))"%src_sele)==1:
                    # position hydro based on location of the carbonyl
                    angle = cmd.get_dihedral("(%s & name C)"%frag_name,
                                             "(%s & name CA)"%frag_name,
                                             "(%s & name N)"%frag_name,
                                             tmp_sele1)
                    cmd.set_dihedral("(%s & name C)"%frag_name,
                                     "(%s & name CA)"%frag_name,
                                     "(%s & name N)"%frag_name,
                                     "(%s & name H)"%frag_name,180.0+angle)
                    cmd.delete(tmp_sele1)

            # add c-cap (if appropriate)
            if self.c_cap in [ 'amin', 'nmet' ]:
                if not cmd.count_atoms("elem N & !(%s) & (bto. (name C & (%s))) & !resn NME+NHH"%
                                       (src_sele,src_sele)):
                    if cmd.count_atoms("name C & (%s)"%(frag_name))==1:
                        if self.c_cap == 'amin':
                            editor.attach_amino_acid("name C & (%s)"%(frag_name), 'nhh', _self=cmd)
                        elif self.c_cap == 'nmet':
                            editor.attach_amino_acid("name C & (%s)"%(frag_name), 'nme', _self=cmd)
                        if cmd.count_atoms("hydro & bound_to (name N & bound_to (name C & (%s)))"%frag_name):
                            cmd.h_fix("name N & bound_to (name C & (%s))"%frag_name)
                        # trim hydrogens
                        if (self.hyd == 'none'):
                            cmd.remove("("+frag_name+" and hydro)")
                        elif (self.hyd == 'auto'):
                            if cmd.count_atoms("("+src_sele+") and hydro")==0:
                                cmd.remove("("+frag_name+" and hydro)")

            # add n-cap (if appropriate)
            if self.n_cap in [ 'acet' ]:
                if not cmd.count_atoms("elem C & !(%s) & (bto. (name N & (%s))) & !resn ACE "%
                                       (src_sele,src_sele)):
                    if cmd.count_atoms("name N & (%s)"%(frag_name))==1:
                        if self.n_cap == 'acet':
                            editor.attach_amino_acid("name N & (%s)"%(frag_name), 'ace', _self=cmd)
                        if cmd.count_atoms("hydro & bound_to (name N & bound_to (name C & (%s)))"%frag_name):
                            cmd.h_fix("name N & (%s)"%frag_name)
                        # trim hydrogens
                        if (self.hyd == 'none'):
                            cmd.remove("("+frag_name+" and hydro)")
                        elif (self.hyd == 'auto'):
                            if cmd.count_atoms("("+src_sele+") and hydro")==0:
                                cmd.remove("("+frag_name+" and hydro)")




        cartoon = (cmd.count_atoms("(%s & name CA & rep cartoon)"%src_sele)>0)
        sticks = (cmd.count_atoms("(%s & name CA & rep sticks)"%src_sele)>0)

        cmd.delete(obj_name)
        key = rot_type
        lib = None
        if self.dep == 'dep':
            try:
                result = cmd.phi_psi("%s"%src_sele)
                if len(result)==1:
                    (phi,psi) = list(result.values())[0]
                    (phi,psi) = (int(10*round(phi/10)),int(10*(round(psi/10))))
                    key = (rot_type,phi,psi)
                    if key not in self.dep_library:
                        (phi,psi) = (int(20*round(phi/20)),int(20*(round(psi/20))))
                        key = (rot_type,phi,psi)
                        if key not in self.dep_library:
                            (phi,psi) = (int(60*round(phi/60)),int(60*(round(psi/60))))
                            key = (rot_type,phi,psi)
                    lib = self.dep_library.get(key,None)
            except:
                pass
        if lib is None:
            key = rot_type
            lib = self.ind_library.get(key,None)
            if (lib is not None) and self.dep == 'dep':
                print(' Mutagenesis: no phi/psi, using backbone-independent rotamers.')
        if lib is not None:
            state = 1
            for a in lib:
                cmd.create(obj_name,frag_name,1,state)
                if state == 1:
                    cmd.select(mut_sele,"(byres (%s like %s))"%(obj_name,src_sele))
                if rot_type=='PRO':
                    cmd.unbond("(%s & name N)"%mut_sele,"(%s & name CD)"%mut_sele)
                for b in a.keys():
                    if b!='FREQ':
                        cmd.set_dihedral("(%s & n;%s)"%(mut_sele,b[0]),
                                         "(%s & n;%s)"%(mut_sele,b[1]),
                                         "(%s & n;%s)"%(mut_sele,b[2]),
                                         "(%s & n;%s)"%(mut_sele,b[3]),
                                         a[b],state=state)
                    else:
                        cmd.set_title(obj_name,state,"%1.1f%%"%(a[b]*100))
                if rot_type=='PRO':
                    cmd.bond("(%s & name N)"%mut_sele,"(%s & name CD)"%mut_sele)
                state = state + 1
            cmd.delete(frag_name)
            print(" Mutagenesis: %d rotamers loaded."%len(lib))
            if self.bump_check:
                cmd.delete(bump_name)
                cmd.create(bump_name,
                "(((byobj %s) within 6 of (%s and not name N+C+CA+O+H+HA)) and (not (%s)))|(%s)"%
                           (src_sele,mut_sele,src_sele,mut_sele),singletons=1)
                cmd.color("gray50",bump_name+" and elem C")
                cmd.set("seq_view",0,bump_name,quiet=1)
                cmd.hide("everything",bump_name)
                if ((cmd.select(tmp_sele1, "(name N & (%s in (neighbor %s)))"%
                                (bump_name,src_sele)) == 1) and
                    (cmd.select(tmp_sele2, "(name C & (%s in %s))"%
                                (bump_name,mut_sele)) == 1)):
                    cmd.bond(tmp_sele1,tmp_sele2)
                if ((cmd.select(tmp_sele1,"(name C & (%s in (neighbor %s)))"%
                                (bump_name,src_sele)) == 1) and
                    (cmd.select(tmp_sele2,"(name N & (%s in %s))"%
                                (bump_name,mut_sele)) == 1)):
                    cmd.bond(tmp_sele1,tmp_sele2)
                cmd.delete(tmp_sele1)
                cmd.delete(tmp_sele2)

                cmd.protect("%s and not (%s in (%s and not name N+C+CA+O+H+HA))"%
                            (bump_name,bump_name,mut_sele))
                cmd.sculpt_activate(bump_name)
                cmd.show("cgo",bump_name)
                # draw the bumps
                cmd.set("sculpt_vdw_vis_mode",1,bump_name)
                state = 1
                score_best = 1e6
                for a in lib:
                    score = cmd.sculpt_iterate(bump_name, state, 1)
                    self.bump_scores.append(score)
                    if score < score_best:
                        state_best = state
                        score_best = score
                    state = state + 1
            cmd.delete(mut_sele)
        else:
            cmd.create(obj_name,frag_name,1,1)
            print(" Mutagenesis: no rotamers found in library.")
        cmd.set("seq_view",0,obj_name,quiet=1)
        pymol.util.cbaw(obj_name, _self=cmd)
        cmd.hide("("+obj_name+")")
        cmd.show(self.rep,obj_name)
        cmd.show('lines',obj_name) #neighbor  always show lines
        if cartoon:
            cmd.show("cartoon",obj_name)
        if sticks:
            cmd.show("sticks",obj_name)
        cmd.set('auto_zoom',auto_zoom,quiet=1)
        cmd.delete(frag_name)
        cmd.frame(state_best)

        # this might be redundant if frame(state_best) changed the state
        self.do_state(state_best)

        cmd.unpick()
        cmd.feedback("pop")

    def do_state(self,state):
        cmd=self.cmd
        if cmd.get("sculpting")=="on":
            names = cmd.get_names("all_objects")
            if (bump_name in names) and (obj_name in names):
                cmd.update(bump_name,obj_name)
        if self.bump_scores:
            print(' Rotamer %d/%d, strain=%.2f' % (state,
                    cmd.count_states(obj_name), self.bump_scores[state - 1]))

        # update hbonds (polar contacts)
        cmd.delete(tmp_hbonds)
        cmd.distance(tmp_hbonds,
                "?{}".format(obj_name),
                "(byobj ?{}) and not ?{}".format(src_sele, src_sele),
                mode=2, state1=state, state2=1)

    def do_select(self,selection):
        print("Selected!")
        cmd=self.cmd
        if (obj_name in cmd.get_names()):
            if cmd.count_atoms("(%s) and (%s)"%(obj_name,selection)):
                cmd.deselect()
                return 1
        if self.status!=0:
            cmd.delete(obj_name)
        cmd.select(src_sele,selection)
        cmd.unpick()
        cmd.enable(src_sele)
        self.status = 1
        self.error = None
        self.do_library()
        cmd.delete(selection)
        cmd.refresh_wizard()
        cmd.deselect()
        return 1

    def do_pick(self,bondFlag):
        print("Picked!")
        cmd=self.cmd
        if bondFlag:
            self.error = "Error: please select an atom, not a bond."
            print(self.error)
        else:
            if self.status!=0:
                cmd.delete(obj_name)
            cmd.select(src_sele,"(byres pk1)")
            cmd.unpick()
            cmd.enable(src_sele)
            self.status = 1
            self.error = None
            self.do_library()
        cmd.refresh_wizard()
