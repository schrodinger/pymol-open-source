#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (C) 2008 DeLano Scientific LLC
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

import cmd as cmd_module
from cmd import _cmd, lock, unlock, Shortcut, \
     _feedback, fb_module, fb_mask, \
     DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error, \
     is_list, safe_list_eval, is_string

import string
import traceback

def model_to_sdf_list(model):
    from chempy import io

    sdf_list = io.mol.toList(model)
    fixed = []
    restrained = []
    at_id = 1
    for atom in model.atom:
        if atom.flags & 4:
            if hasattr(atom,'ref_coord'):
                restrained.append( [at_id,atom.ref_coord])
        if atom.flags & 8:
            fixed.append(at_id)
        at_id = at_id + 1
    fit_flag = 1
    if len(fixed):
        fit_flag = 0
        sdf_list.append(">  <FIXED_ATOMS>\n")
        sdf_list.append("+ ATOM\n");
        for ID in fixed:
            sdf_list.append("| %4d\n"%ID)
        sdf_list.append("\n")
    if len(restrained):
        fit_flag = 0
        sdf_list.append(">  <RESTRAINED_ATOMS>\n")
        sdf_list.append("+ ATOM    MIN    MAX F_CONST         X         Y         Z\n")
        for entry in restrained:
            xrd = entry[1]
            sdf_list.append("| %4d %6.3f %6.3f %6.3f %10.4f %10.4f %10.4f\n"%
                            (entry[0],0,0,3,xrd[0],xrd[1],xrd[2]))
        sdf_list.append("\n")
    sdf_list.append("$$$$\n")
#    for line in sdf_list:
#        print line,
    return (fit_flag, sdf_list)

class CleanJob:
    def __init__(self,self_cmd,sele,message=None,state=-1):
        self.cmd = self_cmd
        if message != None:
            self.cmd.wizard('message',message)
        if state<1:
            state = self_cmd.get_state()
        # this code will moved elsewhere
        ok = 1
        try:
            from freemol import mengine
        except:
            ok = 0
            print "Error: Unable to import module freemol.mengine"
        if ok:
            if not mengine.validate():
                ok = 0
                print "Error: Unable to validate freemol.mengine"
        if not ok:
            pass
            # we can't call warn because this is the not the tcl-tk gui thread
            # warn("Please be sure that FreeMOL is correctly installed.")
        else:
            obj_list = self_cmd.get_object_list("bymol ("+sele+")")
            ok = 0
            result = None
            if len(obj_list)==1:
                obj_name = obj_list[0]
                self_cmd.sculpt_deactivate(obj_name) 
                # eliminate all sculpting information for object
                self.cmd.sculpt_purge()
                self.cmd.set("sculpting",0)
                state = self_cmd.get_state()
#                sdf_list = io.mol.toList(self_cmd.get_model(obj_name,state=state)) + ["$$$$\n"]
                if self_cmd.count_atoms(obj_name+" and flag 2"): # any atoms restrained?
                    self_cmd.reference("validate",obj_name,state) # then we have reference coordinates
                (fit_flag, sdf_list) = model_to_sdf_list(self_cmd.get_model(obj_name,state=state))
                result = mengine.run(string.join(sdf_list,''))
                if result != None:
                    if len(result):
                        clean_sdf = result[0]
                        clean_mol = clean_sdf.split("$$$$")[0]
                        if len(clean_mol):
                            clean_name = "builder_clean_tmp"
                            self_cmd.read_molstr(clean_mol, clean_name, zoom=0)
                            # need to insert some error checking here
                            self_cmd.set("retain_order","1",clean_name)
                            if fit_flag:
                                self_cmd.fit(clean_name, obj_name, matchmaker=4,
                                             mobile_state=1, target_state=state)
                            self_cmd.update(obj_name, clean_name, matchmaker=0,
                                            source_state=1, target_state=state)
                            self_cmd.sculpt_activate(obj_name) 
                            self_cmd.sculpt_deactivate(obj_name) 
                            self_cmd.delete(clean_name)
                            ok = 1

            if not ok:
                # we can't call warn because this is the not the tcl-tk gui thread
                print "Cleanup failed.  Invalid input or software malfuction?"
                if result != None:
                    if len(result)>1:
                        print result[1]
        if message!=None:
            self_cmd.set_wizard()

def clean(selection, present='', state=-1, fix='', restrain='', method='mmff', async=0, _self=cmd_module):
    self_cmd = _self

    clean1_sele = "clean1_tmp"
    clean2_sele = "clean2_tmp"
    clean_obj = "clean_obj"
    r = DEFAULT_SUCCESS

    if self_cmd.select(clean1_sele,selection)>0:

        try:
            if present=='':
                self_cmd.select(clean2_sele," byres (byres ("+selection+") extend 1)") # go out 2 residues
            else:
                self_cmd.select(clean2_sele, clean1_sele+" or ("+present+")")

            self_cmd.create(clean_obj, clean2_sele, zoom=0, source_state=state,target_state=1)

            self_cmd.flag(3,clean_obj+" in ("+clean2_sele+" and not "+clean1_sele+")","set") # fix nearby atoms

            self_cmd.h_add(clean_obj) # fill any open valences

            CleanJob(self_cmd, clean_obj, state, message='Cleaning structure.  Please wait...')

            self_cmd.set_wizard()
            self_cmd.update(clean1_sele, clean_obj, 
                            source_state=1, target_state=state)

            self_cmd.delete(clean_obj)
            self_cmd.delete(clean1_sele)
            self_cmd.delete(clean2_sele)
        except:
            traceback.print_exc()
    return r
