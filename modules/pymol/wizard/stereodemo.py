from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types

saved = {}

class Stereodemo(Wizard):

    def launch(self,name,pretty_name=None):
	demo = DemoInfo()
        if self.last:
            if hasattr(demo,self.last):
                getattr(demo,self.last)(cleanup=1)
        if hasattr(demo,name):
	    cmd.delete("all")
	    if pretty_name != None:
		cmd.do("_ wizard message, Please wait while the %s example loads..., dismiss=0"%pretty_name)
            self.message = demo.message_dict.get(name,None)
	    cmd.refresh_wizard()
            self.last = name
            demo_fn = getattr(demo,name)
            t = threading.Thread(target=demo_fn)
            t.setDaemon(1)
            t.start()
        else:
            self.last = None
	saved['last']=self.last

    def __init__(self,*arg,**kw):
        self.message = []
        self.last = None
        cmd.set("use_display_lists","on")
        cmd.full_screen("off")
        if not  ("mono" in kw.keys()):
            cmd.stereo("on")
        cmd.set("sphere_mode","5")
        if saved.has_key('last'):
            self.last = saved['last']
        if len(arg):
            self.launch(arg[0])
        else:
            self.launch("cartoon")

    def get_prompt(self):
        saved['last']=self.last
        self.prompt = self.message
        return self.prompt

    def get_panel(self):
        return [
            [ 1, 'Structural Biology', '' ],
            [ 2, 'X-ray Crystallography', 'cmd.get_wizard().launch("roving_density")'],
            [ 2, 'Electron Tomography', 
	      'cmd.get_wizard().launch("electomo", "Electron Tomography")'],
            [ 1, 'Drug Discovery', '' ],	    
            [ 2, 'Medicinal Chemistry', 
	      'cmd.get_wizard().launch("medchem","Medicinal Chemistry")'],
            [ 2, 'Computational Chemistry', 'cmd.get_wizard().launch("electro","Computational Chemistry")'],
            [ 1, 'Presentation Graphics', '' ],	    
            [ 2, 'Molecular Animation', 'cmd.get_wizard().launch("animate","Molecular Animation")'],         
            [ 2, 'Multiprocessor Raytracing', 'cmd.get_wizard().launch("ray")'],
	    [ 1, 'Bioinformatics', ''],
            [ 2, 'Structure Alignments', 'cmd.get_wizard().launch("structure","Structure Alignment")'],
            [ 2, 'Homology Modeling',
	      'cmd.get_wizard().launch("homology","Homology Modeling")'],
	    [ 1, 'Science Education', ''],
            [ 2, 'Interactive Modeling', 'cmd.get_wizard().launch("sculpt")'],
            [ 1, 'Configuration', ''],
	    [ 2, 'Toggle Fullscreen', 
	      'cmd.full_screen(apply(lambda x:{ "off":"on", "on":"off"}[x],(cmd.get("full_screen"),)))'],
	    [ 2, 'Toggle Stereo 3D', 
	      'cmd.stereo(apply(lambda x:{ "off":"on", "on":"off"}[x],(cmd.get("stereo"),)))'],
#            [ 2, 'End Demonstration', 'cmd.get_wizard().launch("finish")' ],
#            [ 2, 'Swap Left/Right Stereo', 'stereo swap'],
            ]

from pymol.vfont import plain
from pymol.cgo import *
import string
import traceback
from pymol import util
import threading

class DemoInfo:

    message_dict = {
        'roving' : [ 
        "Middle-click-and-drag to move...",],
        'roving_density' : [
        "Middle-click-and-drag to move..."],
        'elec' : [
        "CTRL-Middle-Click on color bar to change levels...",],
        'sculpt' : [
        "Control-left-click-and-drag on atom centers to drag atoms...",],
        }

    def get_sess(self,file):
        from chempy import io
	try:
	    file = cmd.exp_path(file)
	    sess = io.pkl.fromFile(file)
	    del sess['wizard']
	    del sess['main']  
	    return sess
	except:
	    traceback.print_exc()
	return None

    def homology(self,cleanup=0):
	if not cleanup:
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/homology.pse"))
	    cmd.do("replace_wizard toggle, Homology Modeling")

    def structure(self,cleanup=0):
	if not cleanup:
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/structure.pse"))
	    cmd.do("replace_wizard toggle, Structure Alignment")
	    cmd.set("seq_view_label_mode",1)
	    cmd.set("seq_view",1)
	else:
	    cmd.set("seq_view_label_mode",0)
	    cmd.set("seq_view",0)

    def medchem(self,cleanup=0):
	if not cleanup:
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/drugdisc.pse"))
	    cmd.set("use_display_lists",0)
	    cmd.do("replace_wizard toggle, Medicinal Chemistry")
	else:
	    cmd.set("sphere_scale",1.0)

    def electomo(self,cleanup=0):
	if not cleanup:
	    cmd.feedback("disable","objectsurface","actions")
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/flagellar.pse"))
	    cmd.set("use_display_lists",0)
	    cmd.set("sweep_mode",3)
	    cmd.set("sweep_angle",3)
	    cmd.rock()
	    cmd.do("replace_wizard toggle, Electron Tomography")
	else:
	    cmd.mstop()
	    cmd.rock(0)

    def electro(self,cleanup=0):
	if not cleanup:
	    cmd.set("use_display_lists",0)
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/electro.pse"))
	    cmd.do("replace_wizard toggle, Computational Chemistry (Electrostatics)")
	    cmd.rock(1)
	else:
	    cmd.rock(0)

    def animate(self,cleanup=0):
	if not cleanup:
	    cmd.set("security",0)
	    cmd.set_session(self.get_sess("$PYMOL_DATA/big_demo/animate.pse"))
	    cmd.rock(1)
	    cmd.set("use_display_lists","on")
	    cmd.set("field_of_view",23)
	    cmd.rock(1)
	    cmd.set("sweep_mode",3)
	    cmd.set("sweep_angle",10)
	    cmd.set("sphere_mode",5)
	    cmd.do("replace_wizard toggle, Molecular Animation")
	else:
	    cmd.set("mesh_width",1)
	    cmd.set("field_of_view",20)
	    cmd.rock(0)
	    cmd.mset()
	    cmd.mstop()

    def rep_old(self,cleanup=0):
        if not cleanup:
            try:
                cmd.set("suspend_updates",1,quiet=1)
                cmd.disable()
                cmd.delete("pept")
                cmd.delete("pept_dist")
                cmd.load("$PYMOL_DATA/demo/pept.pdb")
                cmd.show("sticks","(pept and not i;5:7)")
                cmd.show("surface","(pept and i;5,6)")
                cmd.show("mesh","(pept and i;1,11,12,13)")
                cmd.show("spheres","(pept and i;2,12,9,4 and not n;c,n,o,ca)")
                cmd.show("dots","(i;8)")
                cmd.dist("pept_dist","(pept and i;1&n;OD2)","(pept and i;13&n;OG1)")
                cmd.set("dot_width","2");
            finally:
                cmd.set("suspend_updates",0,quiet=1)
        else:
            cmd.delete("pept")
            cmd.delete("pept_dist")

    def reps(self,cleanup=0):
        rep_list = [ "lines","sticks","spheres","surface","mesh","dots","ribbon","cartoon" ]
        try:
            if not cleanup:
                cmd.disable()
                cmd.set("suspend_updates",1,quiet=1)
                cmd.load("$PYMOL_DATA/demo/pept.pdb","rep1")
                cmd.alter("rep1///1-5+8-13/","ss='S'")
                cmd.cartoon("auto")
                cmd.hide("everything","rep1")
                for a in range(2,9):
                    cmd.create("rep%d"%a,"rep1")
                map(lambda x,y:cmd.show(x,"rep%d"%y),
                     rep_list,
                     range(1,9))
                cmd.reset()
                cmd.zoom("rep1",24)
                util.cbay("rep2")
                util.cbac("rep3")
                util.cbas("rep4")
                util.cbab("rep5")
                util.cbaw("rep6")            
                util.cbay("rep8")


                cmd.set("suspend_updates",0,quiet=1)
                scale=0.5
                for b in range(1,20):
                    cmd.set("suspend_updates",0,quiet=1)
                    cmd.refresh()
                    cmd.set("suspend_updates",1,quiet=1)
                    xt=-3.2
                    yt=1.6
                    for a in range(1,5):
                        cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                        xt=xt+2
                    yt=-yt
                    xt=-3.2
                    for a in range(5,9):
                        cmd.translate([xt*scale,yt*scale,0],object="rep%d"%a,camera=0)
                        xt=xt+2
                for a in range(1,9):
                    cmd.origin("rep%d"%a,object="rep%d"%a)
                cmd.mset("1")
                st = string.join(map(lambda x,y:"rotate angle=-3,object=rep%d,axis=%s;"%(x,y),range(1,9),
                                            ['x','y','x','y','x','y','x','y']))
                cmd.mdo(1,st)
                cmd.set("suspend_updates",0,quiet=1)
                cmd.mplay()

                cgo = []
                axes = [[4.5,0.0,0.0],[0.0,3.0,0.0],[0.0,0.0,3.0]]

                c = 1
                for a in rep_list:
                    ext = cmd.get_extent("rep%d"%c)
                    pos = [(ext[0][0]+ext[1][0])/2,
                             (ext[0][1]+ext[1][1])/2+14,
                             (ext[0][2]+ext[1][2])/2]
                    c = c + 1
                    pos[0]=pos[0]-(measure_text(plain,a,axes)/2)
                    wire_text(cgo,plain,pos,a,axes)
                cmd.set("cgo_line_width",1.5)
                cmd.set("auto_zoom",0)
                cmd.load_cgo(cgo,'reps')
                cmd.set("auto_zoom",1)
            else:
                cmd.delete("rep*")
                cmd.mset()
                cmd.mstop()
        except:
            traceback.print_exc()
            
 
    def roving(self,cleanup=0):
        if not cleanup:
            cmd.load("$PYMOL_DATA/demo/il2.pdb")
            cmd.remove("hydro")
            cmd.disable()
            cmd.enable("il2")
            cmd.set("ribbon_color","blue","il2")
            cmd.set("roving_detail",1)
            cmd.set("roving_origin",1)
            cmd.set("stick_radius",0.12,"il2")
#         cmd.zoom("/il2///16/O")
#         cmd.zoom("center",12)
            cmd.set_view ((\
      0.132852688,   -0.729740858,    0.670686543,\
      -0.228543565,    0.635894477,    0.737154961,\
      -0.964425683,   -0.251212329,   -0.082298420,\
      0.000062190,    0.000183226,  -58.861488342,\
      13.349151611,   -1.565427899,   22.383148193,\
      55.259441376,   63.259449005,    0.000000000 ))
        else:
            cmd.delete("il2")
            cmd.set("roving_detail",0)
            cmd.refresh()
            cmd.delete("rov_*")
            
    def roving_density(self,cleanup=0):
        if not cleanup:
            try:
                cmd.load("$PYMOL_DATA/demo/il2.pdb")
                cmd.set("suspend_updates",1,quiet=1)
                cmd.remove("hydro")
                cmd.disable()
                cmd.enable("il2")
                cmd.map_new("map","gaussian","0.75","il2")
		cmd.feedback("disable","objectmesh","actions")
                cmd.set("ribbon_color","purple","il2")
                cmd.set("roving_detail",1)
                cmd.set("roving_origin",1)
                cmd.set("stick_radius",0.12,"il2")
                cmd.set("roving_sticks",0)
                cmd.set("roving_polar_contacts",0)
                cmd.set("line_width","3")
                cmd.set("roving_map1_name","map")
                cmd.isomesh("rov_m1","map",9999.0,"il2")
                cmd.color("density","rov_m1")
                
                cmd.set_view ((\
          0.132852688,   -0.729740858,    0.670686543,\
          -0.228543565,    0.635894477,    0.737154961,\
          -0.964425683,   -0.251212329,   -0.082298420,\
          0.000062190,    0.000183226,  -58.861488342,\
          13.349151611,   -1.565427899,   22.383148193,\
          55.259441376,   63.259449005,    0.000000000 ))
            finally:
                cmd.set("suspend_updates",0,quiet=1)
            cmd.refresh()
        else:
            cmd.set("roving_detail",0)
            cmd.set("roving_map1_name","")
            cmd.set("roving_polar_contacts",7)
            cmd.set("roving_sticks",6)
            cmd.delete("il2")
            cmd.delete("map")
            cmd.set("line_width",1.5)
            cmd.refresh()
            cmd.set("roving_detail",0)
            cmd.delete("rov_*")
            cmd.sync()
            
    def cartoon(self,cleanup=0):
        if not cleanup:
            try:
                cmd.set("suspend_updates",1,quiet=1)
                cmd.disable()
                cmd.delete("1tii")      
                cmd.load("$PYMOL_DATA/demo/1tii.pdb")
                cmd.hide("(1tii)")
                cmd.show("cartoon","1tii")
                cmd.zoom("1tii")
                cmd.spectrum("count","rainbow","1tii////ca")
                cmd.set("cartoon_highlight_color","grey50","1tii")
                cmd.set("cartoon_fancy_helices",1,"1tii")
            finally:
                cmd.set("suspend_updates",0,quiet=1)
            cmd.refresh()
        else:
            cmd.delete("1tii")

    def elec(self,cleanup=0):
        if not cleanup:
            cmd.disable()
            cmd.delete("pept")
            cmd.delete("e_pot")
            cmd.delete("e_lvl")
            cmd.load("$PYMOL_DATA/demo/pept.pkl")
            cmd.hide("(pept)")
            cmd.show("surface","pept")
            cmd.set("coulomb_dielectric",80.0)
            cmd.map_new("e_pot","coulomb",1.0,"pept",5)
            cmd.ramp_new("e_lvl","e_pot",[-3.6,-1.6,0.4])
            cmd.set("surface_color","e_lvl","pept")
            cmd.refresh()
        else:
            cmd.delete("pept")
            cmd.delete("e_pot")
            cmd.delete("e_lvl")
            
    def trans(self,cleanup=0):
        if not cleanup:
            try:
                cmd.set("suspend_updates",1,quiet=1)
                cmd.disable()
                cmd.delete("trans")
                cmd.load("$PYMOL_DATA/demo/pept.pdb","trans")
                cmd.hide("(trans)")
                cmd.show("surface","trans")
                cmd.show("sticks","trans")
                cmd.set("surface_color","white","trans")
                cmd.set("transparency",0.5,"trans")
                cmd.zoom("trans")
            finally:
                cmd.set("suspend_updates",0,quiet=1)
            cmd.refresh()
        else:
            cmd.delete("trans")

    def ray(self,cleanup=0):
        if not cleanup:
            cmd.set("suspend_updates",1,quiet=1)
            cmd.disable()
            cmd.delete("ray")
	    cmd.set("sphere_mode",5)
	    cmd.set("sphere_scale",1.0)
            cmd.load("$PYMOL_DATA/demo/il2.pdb","ray")
            cmd.remove("(ray and hydro)")
            cmd.hide("lines","ray")
            cmd.show("spheres","ray")
            cmd.orient("ray")
            cmd.turn("x",90)
            util.ray_shadows('heavy')
	    cmd.mstop()
	    cmd.rock(0)
            cmd.set("suspend_updates",0,quiet=1)
            cmd.refresh()
            cmd.do("ray")
        else:
            cmd.delete("ray")
            
    def finish(self,cleanup=0):
        cmd.do("_ wizard")

    def sculpt(self,cleanup=0):
        if not cleanup:
            cmd.set("suspend_updates",1,quiet=1)
            cmd.disable()
            cmd.delete("sculpt")
	    cmd.set("sphere_scale","1.0")
	    cmd.set("sphere_mode",5)
            cmd.load("$PYMOL_DATA/demo/pept.pdb","sculpt")
            cmd.hide("lines","sculpt")
#            cmd.show("sticks","sculpt")
            cmd.show("spheres","sculpt")
#            cmd.set("sphere_transparency","0.75","sculpt")
#            cmd.set("sphere_color","grey","sculpt")
            cmd.frame(1)
            cmd.set("auto_sculpt",1)
            cmd.set("sculpting",1)
            cmd.sculpt_activate("sculpt")
	    cmd.set("sculpting_cycles","100")
            cmd.do("edit_mode")
            cmd.set("valence","0.05")
            cmd.set("suspend_updates",0,quiet=0)
	    cmd.sculpt_iterate("sculpt")
	    cmd.alter_state(1,"sculpt","x=x*1.5;y=y*0.1;z=z*1.5")
	    cmd.zoom()

            cmd.unpick()
        else:
            cmd.set("valence","0")
            cmd.set("sculpting",0)
            cmd.set("auto_sculpt",0)
            cmd.delete("sculpt")
            cmd.mouse()


if __name__=='pymol':
   pymol.wizard.stereodemo = Stereodemo()
   cmd.set_wizard(pymol.wizard.stereodemo)
   cmd.stereo("on")
   cmd.set("max_threads",4)
