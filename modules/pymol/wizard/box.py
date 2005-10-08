# This wizard contributed by Ezequiel "Zac" Panepucci 011114
# modified by Warren L. DeLano

from pymol.wizard import Wizard
from pymol import cmd
import pymol
from chempy.models import Indexed
from chempy import Atom,Bond
from chempy.cpv import add, sub, cross_product, scale, dot_product
from chempy.cpv import normalize, project, remove_component, reverse
from pymol.cgo import *

pseudo_atoms = [
   ['X', 'PLN1' , '1', 'PSDO', [0.0,0.0,0.0]],
   ['X', 'PLN2' , '1', 'PSDO', [10.0,0.0,0.0]],
   ['X', 'PLN3' , '1', 'PSDO', [10.0,10.0,0.0]],
   ['X', 'PLN4' , '1', 'PSDO', [10.0,10.0,4]]]

default_mode = 'box'

class Box(Wizard):

    atom=None
    messages=1
    labeling=1
    obj_name=None


    def __init__(self):

        Wizard.__init__(self)
        
        model = Indexed()

        # append the atoms onto it 

        origin = cmd.get_view()[12:15]
        
        for a in pseudo_atoms:
            new_atom = Atom()
            (new_atom.symbol, new_atom.name, new_atom.resi, new_atom.resn, new_atom.coord) = a
            new_atom.coord[0] = new_atom.coord[0] + origin[0]
            new_atom.coord[1] = new_atom.coord[1] + origin[1]
            new_atom.coord[2] = new_atom.coord[2] + origin[2]
            model.atom.append(new_atom)

        cmd.delete("box_points") # just in case
        cmd.load_model(model,"box_points")

        cmd.color("green","box_points`1")
        cmd.color("green","box_points`2")
        cmd.color("red","box_points`3")
        cmd.color("blue","box_points`4")
        cmd.as("nb_spheres","box_points")
        cmd.set("nonbonded_size","0.5","box_points")
        self.points = 1

        self.mode = default_mode
        self.modes = [
            'box',
            'walls', 
            'plane',
            ]

        self.mode_name = {
            'box':'Box',
            'walls':'Walls',
            'plane':'Plane',
            }

        smm = []
        smm.append([ 2, 'Box Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm


        cmd.unpick()

    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.status = 0
        cmd.refresh_wizard()
        
    def get_prompt(self):
        self.prompt = []
        return self.prompt

    def toggle_points(self):
        self.points = not self.points
        if(self.points):
            cmd.enable("box_points")
        else:
            cmd.disable("box_points")

    def update_plane(self):

        model = cmd.get_model("box_points")

        p = model.atom[0].coord
        
        d10 = sub(model.atom[1].coord, p)
        d20 = sub(model.atom[2].coord, p)
        d30 = sub(model.atom[3].coord, p)

        x10_20 = cross_product(d10,d20)
        if dot_product(d30,x10_20)<0.0:
            p = model.atom[1].coord
            d10 = sub(model.atom[0].coord, p)
            d20 = sub(model.atom[2].coord, p)
            d30 = sub(model.atom[3].coord, p)

        n10_20 = normalize(x10_20)
        n10 = normalize(d10)

        d100 = d10
        d010 = remove_component(d20, n10)
        d001 = project(d30, n10_20)

        n100 = normalize(d100)
        n010 = normalize(d010)
        n001 = normalize(d001)

        f100 = reverse(n100)
        f010 = reverse(n010)
        f001 = reverse(n001)
        
        p000 = p
        p100 = add(p,d100)
        p010 = add(p,d010)
        p001 = add(p,d001)
        p110 = add(p100, d010)
        p011 = add(p010, d001)
        p101 = add(p100, d001)
        p111 = add(p110, d001)
                
        obj = []

        if self.mode == 'box': # standard box

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(f001)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p010)
            obj.append(VERTEX); obj.extend(p100)
            obj.append(VERTEX); obj.extend(p110)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n001)
            obj.append(VERTEX); obj.extend(p001)
            obj.append(VERTEX); obj.extend(p101)
            obj.append(VERTEX); obj.extend(p011)
            obj.append(VERTEX); obj.extend(p111)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(f010)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p100)
            obj.append(VERTEX); obj.extend(p001)
            obj.append(VERTEX); obj.extend(p101)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n010)
            obj.append(VERTEX); obj.extend(p010)
            obj.append(VERTEX); obj.extend(p011)
            obj.append(VERTEX); obj.extend(p110)
            obj.append(VERTEX); obj.extend(p111)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(f100)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p001)
            obj.append(VERTEX); obj.extend(p010)
            obj.append(VERTEX); obj.extend(p011)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n100)
            obj.append(VERTEX); obj.extend(p100)
            obj.append(VERTEX); obj.extend(p110)
            obj.append(VERTEX); obj.extend(p101)
            obj.append(VERTEX); obj.extend(p111)
            obj.append(END)

            model.atom[0].coord = p000
            model.atom[1].coord = p100
            model.atom[2].coord = add(p010, scale(d100,0.5))
            model.atom[3].coord = add(add(p001, scale(d010,0.5)),d100)

        elif self.mode=='walls':

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n001)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p100)
            obj.append(VERTEX); obj.extend(p010)
            obj.append(VERTEX); obj.extend(p110)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n010)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p001)
            obj.append(VERTEX); obj.extend(p100)
            obj.append(VERTEX); obj.extend(p101)
            obj.append(END)

            obj.extend([ BEGIN, TRIANGLE_STRIP ])
            obj.append(NORMAL); obj.extend(n100)
            obj.append(VERTEX); obj.extend(p000)
            obj.append(VERTEX); obj.extend(p010)
            obj.append(VERTEX); obj.extend(p001)
            obj.append(VERTEX); obj.extend(p011)
            obj.append(END)

            model.atom[0].coord = p000
            model.atom[1].coord = p100
            model.atom[2].coord = p010
            model.atom[3].coord = p001

        cmd.load_model(model, '_tmp', zoom=0)
        cmd.update("box_points","_tmp")
        cmd.delete("_tmp")
        
        
        # then we load it into PyMOL

        cmd.delete("box_cgo")
        cmd.load_cgo(obj,'box_cgo',zoom=0)
        
    def get_panel(self):
        
        return [
            [ 1, 'Box Wizard',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 2, 'Toggle Points','cmd.get_wizard().toggle_points()'],
            [ 2, 'Update','cmd.get_wizard().update_plane()'],
            
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def do_pick(self,bondFlag):
        pass

