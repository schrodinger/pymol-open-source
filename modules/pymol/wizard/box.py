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
import math

pseudo_atoms = [
   ['X', 'PLN1' , '1', 'PSDO', [0.0,0.0,0.0]],
   ['X', 'PLN2' , '1', 'PSDO', [10.0,0.0,0.0]],
   ['X', 'PLN3' , '1', 'PSDO', [10.0,10.0,0.0]],
   ['X', 'PLN4' , '1', 'PSDO', [10.0,10.0,4]]]

default_mode = 'box'

default_name = 'box'

class Box(Wizard):

    atom=None
    messages=1
    labeling=1
    obj_name=None


    def __init__(self,_self=cmd):
        Wizard.__init__(self,_self)

        self.editing_name = 0
        self.copying = 0

        self.points_name = ''
        self.set_name(default_name)

        self.mode = default_mode
        self.modes = [
            'box',
            'walls',
            'plane',
	    'quad',
            ]

        self.mode_name = {
            'box':'Box',
            'walls':'Walls',
            'plane':'Plane',
	    'quad':'Quad',
            }

        smm = []
        smm.append([ 2, 'Box Mode', '' ])
        for a in self.modes:
            smm.append([ 1, self.mode_name[a], 'cmd.get_wizard().set_mode("'+a+'")'])
        self.menu['mode']=smm

        self.update_box()

    def set_mode(self,mode):
        if mode in self.modes:
            self.mode = mode
        self.status = 0
        self.update_box()
        self.cmd.refresh_wizard()

    def get_prompt(self):
        self.prompt = []
        return self.prompt

    def toggle_points(self):
        if self.points_name in self.cmd.get_names(enabled_only=1):
            self.cmd.disable(self.points_name)
        else:
            self.cmd.enable(self.points_name)

    def delete_box(self):
        self.cmd.delete(self.point_name)
        self.cmd.delete(self.cgo_name)

    def edit_name(self,copying=0):
        self.editing_name = 1
        self.copying = copying
        self.new_name = ''
        self.cmd.refresh_wizard()

    def auto_position(self,fract=0.75,size=1.0):

        if self.points_name in self.cmd.get_names():

            if fract<0.5: fract = 1.0 - fract
            if fract>0.995: fract = 0.995 # minimum 1% distance from clipping planes

            fov = self.cmd.get_setting_float("field_of_view")

            tan_half_fov = math.tan(math.pi*fov/360.0)

            model = self.cmd.get_model(self.points_name)

            view = self.cmd.get_view()

            # locate box 1/4 and 3/4 distance from clipping plane

            one_minus_fract = 1.0 - fract

            plane_z1 = - (view[11] + (one_minus_fract*view[15] + fract*view[16]))
            plane_z2 = - (view[11] + (fract*view[15] + one_minus_fract*view[16]))
            normal = [ 0.0, 0.0, 1.0 ]

            # choose the size of the plane (larger than the fov, but not by much)

            plane_size = size*(one_minus_fract*view[15] + fract*view[16])*tan_half_fov

            obj = []
            plane = [
                [ -plane_size, -plane_size, plane_z1 ],
                [ -plane_size,  plane_size, plane_z1 ],
                [  plane_size, -plane_size, plane_z1 ],
                [ -plane_size, -plane_size, plane_z2 ],
                ]

            # then transform plane coordinates into model space

            plane = list(map( lambda p,v=view: [
               v[0] * p[0] + v[1] * p[1] + v[2]* p[2],
               v[3] * p[0] + v[4] * p[1] + v[5]* p[2],
               v[6] * p[0] + v[7] * p[1] + v[8]* p[2]
               ], plane ))

            normal = ( lambda p,v=view:[
               v[0] * p[0] + v[1] * p[1] + v[2]* p[2],
               v[3] * p[0] + v[4] * p[1] + v[5]* p[2],
               v[6] * p[0] + v[7] * p[1] + v[8]* p[2]
               ])(*(normal,))

            plane = list(map( lambda p,v=view: [
               v[12] + p[0], v[13] + p[1], v[14] + p[2]
               ], plane ))

            model.atom[0].coord = plane[0]
            model.atom[1].coord = plane[1]
            model.atom[2].coord = plane[2]
            model.atom[3].coord = plane[3]

            self.cmd.load_model(model, '_tmp', zoom=0)
            self.cmd.update(self.points_name,"_tmp")
            self.cmd.delete("_tmp")

    def set_name(self,name):

        hidden_name = None

        if self.points_name != '':
            if self.points_name in self.cmd.get_names("all"):
                hidden_name = "_"+self.cgo_name
                self.cmd.disable(self.points_name)
                self.cmd.set_name(self.points_name, hidden_name) # hide

        self.name = name
        self.points_name = self.name + "_points"
        self.cgo_name = self.name
        if self.copying and hidden_name is not None:
            self.cmd.copy(self.points_name, hidden_name, zoom=0)
            print("copy")
        else:
            hidden_name = "_"+self.cgo_name
            if hidden_name in self.cmd.get_names("all"):
                self.cmd.set_name(hidden_name, self.points_name)
        self.copying = 0

        if self.points_name not in self.cmd.get_names():
            model = Indexed()
            origin = self.cmd.get_view()[12:15]
            for a in pseudo_atoms:
                new_atom = Atom()
                (new_atom.symbol, new_atom.name, new_atom.resi, new_atom.resn, new_atom.coord) = a
                new_atom.coord[0] = new_atom.coord[0] + origin[0]
                new_atom.coord[1] = new_atom.coord[1] + origin[1]
                new_atom.coord[2] = new_atom.coord[2] + origin[2]
                new_atom.flags = 0x2200000 # set surface ignore flag
                model.atom.append(new_atom)

            self.cmd.load_model(model,self.points_name,zoom=0)
            self.cmd.set("surface_mode",0,self.points_name) # make sure no surface is shown
            self.coord = None
            self.cmd.color("green","%s`1"%self.points_name)
            self.cmd.color("green","%s`2"%self.points_name)
            self.cmd.color("red"  ,"%s`3"%self.points_name)
            self.cmd.color("blue" ,"%s`4"%self.points_name)
            self.cmd.show_as("nb_spheres",self.points_name)

            self.auto_position(0.75,0.5)

        self.cmd.enable(self.points_name)
        self.points_enabled = 1

    def update_box(self):

        if self.points_name in self.cmd.get_names():

            model = self.cmd.get_model(self.points_name)

            self.coord = (
                model.atom[0].coord,
                model.atom[1].coord,
                model.atom[2].coord,
                model.atom[3].coord,
                )

            p = self.coord[0]

            d10 = sub(self.coord[1], p)
            d20 = sub(self.coord[2], p)
            d30 = sub(self.coord[3], p)

            x10_20 = cross_product(d10,d20)
            if self.mode != 'quad':
                if dot_product(d30,x10_20)<0.0:
                    p = model.atom[1].coord
                    d10 = sub(self.coord[0], p)
                    d20 = sub(self.coord[2], p)
                    d30 = sub(self.coord[3], p)

            n10_20 = normalize(x10_20)
            n10 = normalize(d10)

            d100 = d10
            d010 = remove_component(d20, n10)
            if self.mode != 'quad':
                d001 = project(d30, n10_20)
            else:
                d001 = n10_20

            n100 = normalize(d100)
            n010 = normalize(d010)
            n001 = normalize(d001)

            f100 = reverse(n100)
            f010 = reverse(n010)
            f001 = reverse(n001)

            if self.mode == 'quad':
                p000 = p
                p100 = add(p, remove_component(d10,n001))
                p010 = add(p, remove_component(d20,n001))
                p001 = add(p, remove_component(d30,n001))
            else:
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
            elif self.mode=='plane':
                obj.extend([ BEGIN, TRIANGLE_STRIP ])
                obj.append(NORMAL); obj.extend(n001)
                obj.append(VERTEX); obj.extend(p000)
                obj.append(VERTEX); obj.extend(p100)
                obj.append(VERTEX); obj.extend(p010)
                obj.append(VERTEX); obj.extend(p110)
                obj.append(END)
                model.atom[0].coord = p000
                model.atom[1].coord = p100
                model.atom[2].coord = p010
                model.atom[3].coord = add(add(p001, scale(d010,0.5)),scale(d100,0.5))
            elif self.mode=='quad':
                obj.extend([ BEGIN, TRIANGLE_STRIP ])
                obj.append(NORMAL); obj.extend(n001)
                obj.append(VERTEX); obj.extend(p000)
                obj.append(VERTEX); obj.extend(p100)
                obj.append(VERTEX); obj.extend(p010)
                obj.append(VERTEX); obj.extend(p001)
                obj.append(END)
                model.atom[0].coord = p000
                model.atom[1].coord = p100
                model.atom[2].coord = p010
                model.atom[3].coord = p001

            self.cmd.load_model(model, '_tmp', zoom=0)
            self.cmd.update(self.points_name,"_tmp")
            self.cmd.delete("_tmp")

            # then we load it into PyMOL

            self.cmd.delete(self.cgo_name)
            self.cmd.load_cgo(obj,self.cgo_name,zoom=0)
            self.cmd.order(self.cgo_name+" "+self.points_name,sort=1,location='bottom')
            self.cmd.set("nonbonded_size",math.sqrt(dot_product(d10,d10))/10,self.points_name)

    def get_panel(self):

        return [
            [ 1, 'Box Wizard',''],
            [ 3, self.mode_name[self.mode],'mode'],
            [ 2, 'Change Name','cmd.get_wizard().edit_name()'],
            [ 2, 'Copy Box','cmd.get_wizard().edit_name(copying=1)'],
            [ 2, 'Toggle Points','cmd.get_wizard().toggle_points()'],
#            [ 2, 'Update','cmd.get_wizard().update_box()'],
            [ 2, 'Auto-Position (50%)','cmd.get_wizard().auto_position(0.75)'],
            [ 2, 'Auto-Position (99%)','cmd.get_wizard().auto_position(0.99)'],
            [ 2, 'Done','cmd.set_wizard()'],
            ]

    def get_event_mask(self):
        if self.editing_name:
            return Wizard.event_mask_pick + Wizard.event_mask_select + \
                   Wizard.event_mask_scene + Wizard.event_mask_key
        else:
            return Wizard.event_mask_pick + Wizard.event_mask_select + Wizard.event_mask_scene

    def do_scene(self):
        if self.points_name in self.cmd.get_names("objects"):
            if self.coord is None:
                self.update_box()
            else:
                model = self.cmd.get_model(self.points_name)
                coord = (
                    model.atom[0].coord,
                    model.atom[1].coord,
                    model.atom[2].coord,
                    model.atom[3].coord,
                    )
                if self.coord != coord:
                    self.update_box()

    def do_pick(self,bondFlag):
        pass

    def do_key(self,k,x,y,m):
        if k in [8,127]:
            self.new_name = self.new_name[:-1]
        elif k>32:
            self.new_name = self.new_name + chr(k)
        elif k==10 or k==13:
            self.editing_name = 0
            self.new_name = self.new_name.strip()
            if self.new_name == '':
                self.new_name = 'box'
            else:
                self.new_name = self.new_name
            self.set_name(self.new_name)
            self.update_box()
        self.cmd.refresh_wizard()
        return 1

    def get_prompt(self):
        if self.editing_name:
            self.prompt = [ "Enter box name: " + self.new_name ]
        else:
            self.prompt = None
        return self.prompt
