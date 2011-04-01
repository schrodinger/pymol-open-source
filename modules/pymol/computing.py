#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC
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
import threading
import os

def model_to_sdf_list(self_cmd,model):
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
    electro_mode = int(self_cmd.get('clean_electro_mode'))
    if electro_mode == 0:
        fit_flag = 0
        sdf_list.append(">  <ELECTROSTATICS>\n")
        sdf_list.append("+ TREATMENT\n")
        sdf_list.append("| NONE\n")
        sdf_list.append("\n")
    sdf_list.append("$$$$\n")
#    for line in sdf_list:
#        print line,
    return (fit_flag, sdf_list)

def get_energy_from_rec(rec):
    # we really need to replace this with a proper SD parser...
    result = 9999.00
    try:
        rec_list = rec.splitlines()
        read_energy = 0
        for line in rec_list:
            if read_energy == 1:
                result = float(line.strip())
                break
            if line.strip() == '> <MMFF94 energy>':
                read_energy = 1
    except:
        traceback.print_exc()
    return result
    
class CleanJob:
    def __init__(self,self_cmd,sele,state=-1,message=None):
        self.cmd = self_cmd
        if message == '':
            message = None
        if state<1:
            state = self_cmd.get_state()
        # this code will moved elsewhere
        ok = 1
        try:
            from freemol import mengine
        except:
            ok = 0
            print "Error: unable to import freemol.mengine module."
            print "This PyMOL build appears not to include full modeling capabilities."
        if ok:
            if not mengine.validate():
                ok = 0
                print "Error: Unable to validate freemol.mengine"
        if ok:
            if self_cmd.count_atoms(sele) > 999:
                ok = 0
                print "Error: Sorry, clean is currently limited to 999 atoms"
        if not ok:
            pass
            # we can't call warn because this is the not the tcl-tk gui thread
            # warn("Please be sure that FreeMOL is correctly installed.")
        else:
            if message != None:
                self.cmd.do("_ cmd.wizard('message','''%s''')"%message)
            obj_list = self_cmd.get_object_list("bymol ("+sele+")")
            ok = 0
            result = None
            if is_list(obj_list) and (len(obj_list)==1):
                obj_name = obj_list[0]
                self_cmd.sculpt_deactivate(obj_name) 
                # eliminate all sculpting information for object
                self.cmd.sculpt_purge()
                self.cmd.set("sculpting",0)
                state = self_cmd.get_state()
                if self_cmd.count_atoms(obj_name+" and flag 2"): # any atoms restrained?
                    self_cmd.reference("validate",obj_name,state) # then we have reference coordinates
                input_model = self_cmd.get_model(obj_name,state=state)
                (fit_flag, sdf_list) = model_to_sdf_list(self_cmd,input_model)
                input_sdf = string.join(sdf_list,'')
#                print input_sdf
                result = mengine.run(input_sdf)
                if result != None:
                    if len(result):
                        clean_sdf = result[0]
                        clean_rec = clean_sdf.split("$$$$")[0]
                        energy = get_energy_from_rec(clean_rec)
                        if len(clean_rec) and int(energy) != 9999:
                            clean_name = "builder_clean_tmp"
                            self_cmd.set("suspend_updates")
                            try:
                                self_cmd.read_molstr(clean_rec, clean_name, zoom=0)
                                # need to insert some error checking here
                                if clean_name in self_cmd.get_names("objects"):
                                    self_cmd.set("retain_order","1",clean_name)
                                    if fit_flag:
                                        self_cmd.fit(clean_name, obj_name, matchmaker=4,
                                                     mobile_state=1, target_state=state)
                                    self_cmd.push_undo(obj_name)
                                    self_cmd.update(obj_name, clean_name, matchmaker=0,
                                                    source_state=1, target_state=state)
                                    self_cmd.sculpt_activate(obj_name) 
                                    self_cmd.sculpt_deactivate(obj_name)
                                    ok = 1
                            finally:
                                self_cmd.delete(clean_name)
                                self_cmd.unset("suspend_updates")
            if not ok:
                # we can't call warn because this is the not the tcl-tk gui thread
                if result != None:
                    if len(result)>1:
                        print "\n=== mengine errors below === "
                        print result[1].replace("\n\n","\n"),
                        print "=== mengine errors above ===\n"
                failed_file = "cleanup_failed.sdf"
                print "Clean-Error: Structure cleanup failed.  Invalid input or software malfuction?"
                aromatic = 0
                for bond in input_model.bond:
                    if bond.order == 4:
                        aromatic = 1
                try:
                    open(failed_file,'wb').write(input_sdf)
                    print "Clean-Error: Wrote SD file '%s' into the directory:"%failed_file
                    print "Clean-Error: '%s'."%os.getcwd()
                    print "Clean-Error: If you believe PyMOL should be able to handle this structure"
                    print "Clean-Error: then please email that SD file to help@schrodinger.com. Thank you!"
                except IOError:
                    print "Unabled to write '%s"%failed_file
                if aromatic:
                    print "Clean-Warning: Please eliminate aromatic bonds and then try again."                    
        if message!=None:
            self_cmd.do("_ wizard")

def _clean(selection, present='', state=-1, fix='', restrain='',
          method='mmff', async=0, save_undo=1, message=None,
          _self=cmd_module):

    self_cmd = _self

    clean1_sele = "_clean1_tmp"
    clean2_sele = "_clean2_tmp"
    clean_obj = "_clean_obj"
    r = DEFAULT_SUCCESS

    if self_cmd.select(clean1_sele,selection,enable=0)>0:
        try:
            if present=='':
                self_cmd.select(clean2_sele," byres (byres ("+selection+") extend 1)",enable=0) # go out 2 residues
            else:
                self_cmd.select(clean2_sele, clean1_sele+" or ("+present+")",enable=0)

            self_cmd.set("suspend_updates")
            self_cmd.rename(clean2_sele) # ensure identifiers are unique
            self_cmd.create(clean_obj, clean2_sele, zoom=0, source_state=state,target_state=1)
            self_cmd.disable(clean_obj)
            self_cmd.unset("suspend_updates")
            
            self_cmd.flag(3,clean_obj+" in ("+clean2_sele+" and not "+clean1_sele+")","set")
            # fix nearby atoms

            self_cmd.h_add(clean_obj) # fill any open valences

            if message == None:
                at_cnt = self_cmd.count_atoms(clean_obj)
                message = 'Cleaning %d atoms.  Please wait...'%at_cnt

            CleanJob(self_cmd, clean_obj, state, message=message)
            
            self_cmd.push_undo(clean1_sele)
            self_cmd.update(clean1_sele, clean_obj, 
                            source_state=1, target_state=state)

            self_cmd.delete(clean_obj)
            self_cmd.delete(clean1_sele)
            self_cmd.delete(clean2_sele)
        except:
            traceback.print_exc()
    return r

def clean(selection, present='', state=-1, fix='', restrain='',
          method='mmff', async=0, save_undo=1, message=None,
          _self=cmd_module):
    if not int(async):
        return _clean(selection,present,state,fix,restrain,method,async,save_undo,message,_self)
    else:
        try:
            t = threading.Thread(target=_clean,
                             args=(selection,present,state,fix,restrain,
                                   method,async,save_undo,message,_self))
            t.setDaemon(1)
            t.start()
        except:
            traceback.print_exc()
        return 0


