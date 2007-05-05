from pymol.wizard import Wizard
from pymol import cmd
import pymol

class Annotation(Wizard):

    def get_event_mask(self):
        return Wizard.event_mask_scene+Wizard.event_mask_state+Wizard.event_mask_frame
    
    def do_scene(self):
        self.cmd.dirty_wizard()
        
    def do_frame(self,frame):
        self.cmd.dirty_wizard()

    def do_state(self,state):
        self.cmd.dirty_wizard()
            
    def get_prompt(self):
        prompt = []
        if hasattr(pymol.session,'annotation'):
            anno_dict = pymol.session.annotation
            for obj in self.cmd.get_names('objects',1): # enabled objects
                state_dict = anno_dict.get(obj,{})
                state = self.cmd.get_state()
                anno_list = state_dict.get(state,[])
                prompt.extend(anno_list)
        return prompt

    def get_panel(self):
        return [
            [ 1, 'Annotation', '' ],
            [ 2, 'Dismiss', 'cmd.set_wizard()' ]
            ]

import copy
import string
import re

from chempy.sdf import SDF

def load_annotated_sdf(filename, object=None, state=1, discrete=1, _self=cmd):
    pymol=_self._pymol
    cmd=_self
    
    # get object name from file prefix

    if object==None:
        object = re.sub(r"\.[sS][dD][fF]$","",filename)

    # open the SD file

    inp_sdf = SDF(filename)

    # create a persistent place to store the annotations

    if not hasattr(pymol.session,'annotation'):
        pymol.session.annotation = {}

    # create a state-indexed dictionary for this object

    state_dict = {}
    pymol.session.annotation[object] = state_dict

    while 1:

        # get next record 

        sdf_rec = inp_sdf.read()

        # if at end of list, break out of loop
        if not sdf_rec: break

        # get the MOL portion of the record

        mol_list = sdf_rec.get('MOL')

        # load it into PyMOL

        cmd.read_molstr(string.join(mol_list,''),object,
                             state,finish=0,discrete=discrete)

        # populate with tuple containing ordered list of keys
        # and associated data dictionary

        anno_list = [ "\\955"+object ]
        for key in sdf_rec.kees:
            if (key!='MOL'):
                data = sdf_rec.data[key]
                print key,data
                anno_list.append("  \\595%s: \\559%s"%(
                    key,
                    string.join(map(string.strip,sdf_rec.data[key]))))
        state_dict[state] = anno_list

        # increment the state index 

        state = state + 1

    if state > 1:
        cmd.zoom(object)
        cmd.finish_object(object)
