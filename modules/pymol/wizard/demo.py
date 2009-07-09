from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

saved = {}


class Demo(Wizard):

    def launch(self,name):
        return None

    def __init__(self,name=None,_self=cmd):
        Wizard.__init__(self,_self)
        self.message = []
        self.last = None
        if saved.has_key('last'):
            self.last = saved['last']
        if name!=None:
            demo = DemoInfo(_self=_self)
            if self.last:
                if hasattr(demo,self.last):
                    getattr(demo,self.last)(cleanup=1)
            if hasattr(demo,name):
                self.last = name
                demo_fn = getattr(demo,name)
                t = threading.Thread(target=demo_fn)
                t.setDaemon(1)
                t.start()
                self.message = demo.message_dict.get(name,None)
            else:
                self.last = None
            saved['last']=self.last
        
    def get_prompt(self):
        saved['last']=self.last
        self.prompt = self.message
        return self.prompt

    def get_panel(self):
        return [
            [ 1, 'Demonstrations', '' ],
            [ 2, 'Representations', 'replace_wizard demo,reps'],
            [ 2, 'Cartoon Ribbons', 'replace_wizard demo,cartoon'],
            [ 2, 'Roving Detail', 'replace_wizard demo,roving'],         
            [ 2, 'Roving Density', 'replace_wizard demo,roving_density'],
            [ 2, 'Transparency', 'replace_wizard demo,trans'],
            [ 2, 'Ray Tracing', 'replace_wizard demo,ray'],
            [ 2, 'Sculpting', 'replace_wizard demo,sculpt'],
            [ 2, 'Scripted Animation', 'replace_wizard demo,anime'],
            [ 2, 'Electrostatics', 'replace_wizard demo,elec'],
            [ 2, 'CGOs', 'replace_wizard demo,cgo'],
            [ 2, 'Molscript/R3D Input', 'replace_wizard demo,raster3d'],
            [ 2, 'End Demonstration', 'replace_wizard demo,finish' ]
            ]

from pymol.vfont import plain
from pymol.cgo import *
import string
import traceback
from pymol import util
import threading

class DemoInfo:

    def __init__(self,_self=cmd):
        self.cmd=_self
        
    message_dict = {
        'roving' : [ 
        "Middle-Click to rove...         CTRL-SHIFT-Middle-Click to center...",],
        'roving_density' : [
        "Middle-Click to rove...         CTRL-SHIFT-Middle-Click to center...",],
        'elec' : [
        "CTRL-Middle-Click on color bar to change levels...",],
        'sculpt' : [
        "CTRL-Left-Click to drag atoms...       CTRL-Right-Click to rotate bonds...",],
        }
                                          
    def rep_old(self,cleanup=0):
        if not cleanup:
            try:
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.disable()
                self.cmd.delete("pept")
                self.cmd.delete("pept_dist")
                self.cmd.load("$PYMOL_DATA/demo/pept.pdb")
                self.cmd.show("sticks","(pept and not i;5:7)")
                self.cmd.show("surface","(pept and i;5,6)")
                self.cmd.show("mesh","(pept and i;1,11,12,13)")
                self.cmd.show("spheres","(pept and i;2,12,9,4 and not n;c,n,o,ca)")
                self.cmd.show("dots","(i;8)")
                self.cmd.dist("pept_dist","(pept and i;1&n;OD2)","(pept and i;13&n;OG1)")
                self.cmd.set("dot_width","2");
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
        else:
            self.cmd.delete("pept")
            self.cmd.delete("pept_dist")

    def reps(self,cleanup=0):
        rep_list = [ "lines","sticks","spheres","surface","mesh","dots","ribbon","cartoon" ]
        try:
            if not cleanup:
                self.cmd.disable()
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.mset()
                self.cmd.unset("movie_panel")
                self.cmd.load("$PYMOL_DATA/demo/pept.pdb","rep1")
                self.cmd.alter("rep1///1-5+8-13/","ss='S'")
                self.cmd.cartoon("auto")
                self.cmd.hide("everything","rep1")
                for a in range(2,9):
                    self.cmd.create("rep%d"%a,"rep1")
                map(lambda x,y,s=self:s.cmd.show(x,"rep%d"%y),
                     rep_list,
                     range(1,9))
                self.cmd.reset()
                self.cmd.zoom("rep1",24)
                util.cbay("rep2",_self=self.cmd)
                util.cbac("rep3",_self=self.cmd)
                util.cbas("rep4",_self=self.cmd)
                util.cbab("rep5",_self=self.cmd)
                util.cbaw("rep6",_self=self.cmd)            
                util.cbay("rep8",_self=self.cmd)


                self.cmd.set("suspend_updates",0,quiet=1)
                scale=0.5
                for b in range(1,20):
                    self.cmd.set("suspend_updates",0,quiet=1)
                    self.cmd.refresh()
                    self.cmd.set("suspend_updates",1,quiet=1)
                    xt=-3.2
                    yt=1.6
                    for a in range(1,5):
                        self.cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                        xt=xt+2
                    yt=-yt
                    xt=-3.2
                    for a in range(5,9):
                        self.cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                        xt=xt+2
                for a in range(1,9):
                    self.cmd.origin("rep%d"%a,object="rep%d"%a)
                self.cmd.mset("1")
                st = string.join(map(lambda x,y:"rotate angle=-3,object=rep%d,axis=%s;"%(x,y),range(1,9),
                                            ['x','y','x','y','x','y','x','y']))
                self.cmd.mdo(1,st)
                self.cmd.set("suspend_updates",0,quiet=1)
                self.cmd.mplay()

                cgo = []
                axes = [[4.5,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]]

                c = 1
                for a in rep_list:
                    ext = self.cmd.get_extent("rep%d"%c)
                    pos = [(ext[0][0]+ext[1][0])/2,
                             (ext[0][1]+ext[1][1])/2+14,
                             (ext[0][2]+ext[1][2])/2]
                    c = c + 1
                    pos[0]=pos[0]-(measure_text(plain,a,axes)/2)
                    wire_text(cgo,plain,pos,a,axes)
                self.cmd.set("cgo_line_width",1.5)
                self.cmd.set("auto_zoom",0)
                self.cmd.load_cgo(cgo,'reps')
                self.cmd.set("auto_zoom",1)
            else:
                self.cmd.delete("rep*")
                self.cmd.mset()
                self.cmd.mstop()
                self.cmd.set("movie_panel",1)
        except:
            traceback.print_exc()
            
    def raster3d(self,cleanup=0):
        if not cleanup:
            self.cmd.disable()

            self.cmd.set_view( (\
      0.269525230,   -0.492282957,    0.827655137,\
     -0.158114254,   -0.870419860,   -0.466229916,\
      0.949923635,   -0.005200397,   -0.312437057,\
     -0.000086844,    0.000019042, -133.217041016,\
     11.377667427,   21.768899918,    9.270449638,\
    105.029335022,  169.626159668,    0.000000000 ))
            self.cmd.load("$PYMOL_DATA/demo/1hpv.r3d","cgo1")
            self.cmd.zoom("cgo1")
        else:
            self.cmd.delete("cgo1")

    def cgo(self,cleanup=0):
        if not cleanup:
            self.cmd.disable()
            try:
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.do("run $PYMOL_DATA/demo/cgo03.py")
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
        else:
            self.cmd.delete("cgo03")
            self.cmd.mset()
            self.cmd.mstop()
            self.cmd.rewind()

    def anime_old(self,cleanup=0):
        if not cleanup:
            self.cmd.disable()
            self.cmd.delete("arg")
            self.cmd.fragment("arg")
            self.cmd.zoom("arg",2)
            self.cmd.show("sticks","arg")
            self.cmd.feedback('dis','sel','res')
            for a in xrange(1,181):
                try:
                    self.cmd.set("suspend_updates",1,quiet=1)
                    self.cmd.edit("(arg and n;cd)","(arg and n;cg)",quiet=1)
                    self.cmd.torsion("6")
                    self.cmd.unpick()
                    self.cmd.edit("(arg and n;cb)","(arg and n;ca)",quiet=1)
                    self.cmd.torsion("2")
                    self.cmd.unpick()
                finally:
                    self.cmd.set("suspend_updates",0,quiet=1)         
                self.cmd.refresh()
            self.cmd.feedback('ena','sel','res')
        else:
            self.cmd.delete("arg")

    def anime(self,cleanup=0):
        if not cleanup:
            try:
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.disable()
                self.cmd.load("$TUT/1hpv.pdb")
                util.chainbow("1hpv",_self=self.cmd)
                self.cmd.hide("everything","1hpv")
                self.cmd.show("cartoon","1hpv")
                self.cmd.show("sticks","1hpv///200/")
                self.cmd.create("1hpv_a","1hpv//A//")
                self.cmd.set("cartoon_smooth_loops",0,"1hpv_a")
                self.cmd.create("1hpv_b","1hpv//B//")
                self.cmd.set("cartoon_smooth_loops",0,"1hpv_b")
                self.cmd.create("1hpv_l","1hpv///200/")
                util.cbay("1hpv_l",_self=self.cmd)
                self.cmd.delete("1hpv")
                self.cmd.set_view ((\
          0.374249548,   -0.517475128,    0.769516647,\
          -0.214397043,   -0.855623126,   -0.471108317,\
          0.902203023,    0.011330833,   -0.431161582,\
          -0.000023194,   -0.000007302, -125.089942932,\
          11.953758240,   20.323493958,    8.406080246,\
          75.304412842,  189.396347046,    0.000000000 ))
                self.cmd.translate([-20,0,0],object="1hpv_a")
                self.cmd.translate([20,0,0],object="1hpv_b")
                self.cmd.zoom("center",30)
                self.cmd.translate([0,10,00],object="1hpv_l")
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
            self.cmd.refresh()
            for a in range(1,21):
                try:
                    self.cmd.set("suspend_updates",1,quiet=1)
                    self.cmd.translate([1,0,0],object="1hpv_a")
                    self.cmd.translate([-1,0,0],object="1hpv_b")
                    self.cmd.translate([0,-0.5,0],object="1hpv_l")
                finally:
                    self.cmd.set("suspend_updates",0,quiet=1)
                self.cmd.refresh()
            for a in range(1,62):
                self.cmd.turn("y",6)
                self.cmd.move('z',2)
                self.cmd.move('y',-0.12)
                self.cmd.refresh()
                                
        else:
            self.cmd.delete("1hpv_*")

    def roving(self,cleanup=0):
        if not cleanup:
            self.cmd.load("$PYMOL_DATA/demo/il2.pdb")
            self.cmd.remove("hydro")
            self.cmd.disable()
            self.cmd.enable("il2")
            self.cmd.set("ribbon_color","blue","il2")
            self.cmd.set("roving_detail",1)
            self.cmd.set("roving_origin",1)
            self.cmd.set("stick_radius",0.12,"il2")
#         self.cmd.zoom("/il2///16/O")
#         self.cmd.zoom("center",12)
            self.cmd.set_view ((\
      0.132852688,   -0.729740858,    0.670686543,\
      -0.228543565,    0.635894477,    0.737154961,\
      -0.964425683,   -0.251212329,   -0.082298420,\
      0.000062190,    0.000183226,  -58.861488342,\
      13.349151611,   -1.565427899,   22.383148193,\
      55.259441376,   63.259449005,    0.000000000 ))
        else:
            self.cmd.delete("il2")
            self.cmd.set("roving_detail",0)
            self.cmd.refresh()
            self.cmd.delete("rov_*")
            
    def roving_density(self,cleanup=0):
        if not cleanup:
            try:
                self.cmd.load("$PYMOL_DATA/demo/il2.pdb")
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.remove("hydro")
                self.cmd.disable()
                self.cmd.enable("il2")
                self.cmd.map_new("map","gaussian","0.75","il2")
                self.cmd.set("ribbon_color","purple","il2")
                self.cmd.set("roving_detail",1)
                self.cmd.set("roving_origin",1)
                self.cmd.set("stick_radius",0.12,"il2")
                self.cmd.set("roving_sticks",0)
                self.cmd.set("roving_polar_contacts",0)
                self.cmd.set("line_width","3")
                self.cmd.set("roving_map1_name","map")
                self.cmd.isomesh("rov_m1","map",9999.0,"il2")
                self.cmd.color("density","rov_m1")
                
                self.cmd.set_view ((\
          0.132852688,   -0.729740858,    0.670686543,\
          -0.228543565,    0.635894477,    0.737154961,\
          -0.964425683,   -0.251212329,   -0.082298420,\
          0.000062190,    0.000183226,  -58.861488342,\
          13.349151611,   -1.565427899,   22.383148193,\
          55.259441376,   63.259449005,    0.000000000 ))
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
            self.cmd.refresh()
        else:
            self.cmd.set("roving_detail",0)
            self.cmd.set("roving_map1_name","")
            self.cmd.set("roving_polar_contacts",7)
            self.cmd.set("roving_sticks",6)
            self.cmd.delete("il2")
            self.cmd.delete("map")
            self.cmd.set("line_width",1.5)
            self.cmd.refresh()
            self.cmd.set("roving_detail",0)
            self.cmd.delete("rov_*")
            self.cmd.sync()
            
    def cartoon(self,cleanup=0):
        if not cleanup:
            try:
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.disable()
                self.cmd.delete("1tii")      
                self.cmd.load("$PYMOL_DATA/demo/1tii.pdb")
                self.cmd.hide("(1tii)")
                self.cmd.show("cartoon","1tii")
                self.cmd.zoom("1tii")
                self.cmd.spectrum("count","rainbow","1tii////ca")
                self.cmd.set("cartoon_highlight_color","grey50","1tii")
                self.cmd.set("cartoon_fancy_helices",1,"1tii")
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
            self.cmd.refresh()
        else:
            self.cmd.delete("1tii")

    def elec(self,cleanup=0):
        if not cleanup:
            self.cmd.disable()
            self.cmd.delete("pept")
            self.cmd.delete("e_pot")
            self.cmd.delete("e_lvl")
            self.cmd.load("$PYMOL_DATA/demo/pept.pkl")
            self.cmd.hide("(pept)")
            self.cmd.show("surface","pept")
            self.cmd.set("coulomb_dielectric",80.0)
            self.cmd.map_new("e_pot","coulomb",1.0,"pept",5)
            self.cmd.ramp_new("e_lvl","e_pot",[-3.6,-1.6,0.4])
            self.cmd.set("surface_color","e_lvl","pept")
            self.cmd.refresh()
        else:
            self.cmd.delete("pept")
            self.cmd.delete("e_pot")
            self.cmd.delete("e_lvl")
            
    def trans(self,cleanup=0):
        if not cleanup:
            try:
                self.cmd.set("suspend_updates",1,quiet=1)
                self.cmd.disable()
                self.cmd.delete("trans")
                self.cmd.load("$PYMOL_DATA/demo/pept.pdb","trans")
                self.cmd.hide("(trans)")
                self.cmd.show("surface","trans")
                self.cmd.show("sticks","trans")
                self.cmd.set("surface_color","white","trans")
                self.cmd.set("transparency",0.5,"trans")
                self.cmd.zoom("trans")
            finally:
                self.cmd.set("suspend_updates",0,quiet=1)
            self.cmd.refresh()
        else:
            self.cmd.delete("trans")

    def ray(self,cleanup=0):
        if not cleanup:
            self.cmd.set("suspend_updates",1,quiet=1)
            self.cmd.disable()
            self.cmd.delete("ray")
            self.cmd.load("$PYMOL_DATA/demo/il2.pdb","ray")
            self.cmd.remove("(ray and hydro)")
            self.cmd.hide("lines","ray")
            self.cmd.show("spheres","ray")
            self.cmd.orient("ray")
            self.cmd.turn("x",90)
            util.ray_shadows('heavy',_self=self.cmd)
            self.cmd.set("suspend_updates",0,quiet=1)
            self.cmd.refresh()
            self.cmd.do("ray")
        else:
            self.cmd.delete("ray")
            
    def finish(self,cleanup=0):
        self.cmd.do("_ wizard")

    def sculpt(self,cleanup=0):
        if not cleanup:
            self.cmd.set("suspend_updates",1,quiet=1)
            self.cmd.disable()
            self.cmd.delete("sculpt")
            self.cmd.load("$PYMOL_DATA/demo/pept.pdb","sculpt")
            self.cmd.hide("lines","sculpt")
            self.cmd.show("sticks","sculpt")
            self.cmd.show("spheres","sculpt")
            self.cmd.set("sphere_transparency","0.75","sculpt")
            self.cmd.set("sphere_color","grey","sculpt")
            self.cmd.frame(1)
            self.cmd.set("auto_sculpt",1)
            self.cmd.set("sculpting",1)
            self.cmd.sculpt_activate("sculpt")
            self.cmd.do("edit_mode")
            self.cmd.set("valence","0.05")
            self.cmd.set("suspend_updates",0,quiet=0)
            self.cmd.unpick()
        else:
            self.cmd.set("valence","0")
            self.cmd.set("sculpting",0)
            self.cmd.set("auto_sculpt",0)
            self.cmd.delete("sculpt")
            self.cmd.mouse()

