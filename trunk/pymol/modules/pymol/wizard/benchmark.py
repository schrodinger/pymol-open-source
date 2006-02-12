import threading
from pymol.wizard import Wizard
from pymol import cmd
import pymol
import types
import time

def bench_fn(action):
    time.sleep(0.5)
    cmd.do("_ wizard benchmark,%s"%action)

def report(name,value):
    ver = cmd.get_version()[0]
    print "PyMOL %s benchmark: %30s = %10.5f"%(ver,name,value)
    
class Benchmark(Wizard):

    def launch(self,name):
        return None

    def configure(self):
        cmd.reinitialize()
        cmd.set('use_display_lists',1)
        
    def __init__(self,*arg):
        self.gl = 5.0
        self.short_cpu = 8.0
        self.long_cpu = 16.0
        self.message = []
        if len(arg):
            if hasattr(self,arg[0]):
                getattr(self,arg[0])()
                
    def reset(self):
        pass

    def run_all(self):
        self.run_gl()
        self.run_cpu()

    def run_cpu(self):
        self.surface_calculation()
        self.configure()
        self.mesh_calculation()
        self.configure()
        self.ray_trace1()
        self.configure()

    def run_gl(self):
        self.configure()
        self.updates()
        self.configure()
        self.smooth_lines()
        self.configure()
        self.jagged_lines()
        self.configure()
        self.dots()
        self.configure()
        self.sticks()
        self.configure()
        self.surface()
        self.configure()
        self.spheres()
        self.configure()
        self.cartoon()
        self.configure()
        self.blits()
        self.configure()

    def updates(self):
        cmd.fragment("methane")
        cmd.set("antialias",0)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",1)
            cmd.turn("y",1)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('UPDATES_V1',(cnt/elapsed)/100)

    def smooth_lines(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.show("mesh")
        cmd.zoom(complete=1)
        elapsed = 0.0
        cnt = 0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('SMOOTH_LINES_V1',cnt/elapsed)
        
    def jagged_lines(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.show("mesh")
        cmd.set("line_smooth",0)
        cmd.zoom(complete=1)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('JAGGED_LINES_V1',cnt/elapsed)

    def dots(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.hide()
        cmd.show("dots")
        cmd.zoom(complete=1)
        elapsed = 0.0
        cnt = 0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('DOTS_V1',cnt/elapsed)

    def sticks(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.hide()
        cmd.show("sticks")
        cmd.zoom(complete=1)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('STICKS_V1',cnt/elapsed)

    def surface(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.hide()
        cmd.show("surface")
        cmd.zoom(complete=1)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('SURFACE_V1',cnt/elapsed)

    def spheres(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.hide()
        cmd.show("spheres")
        cmd.zoom(complete=1)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('SPHERES_V1',cnt/elapsed)

    def cartoon(self):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.hide()
        cmd.show("cartoon")
        cmd.spectrum("count",selection="name ca")
        cmd.zoom(complete=1)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl: 
            cmd.turn("x",15)
            cmd.turn("y",15)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('CARTOON_V1',cnt/elapsed)

    def blits(self):
        cmd.load("$PYMOL_DATA/demo/pept.pdb")
        cmd.mset("1 x2")
        cmd.set('cache_frames',1)
        cmd.rewind()
        cmd.refresh()
        cmd.turn('x',5)
        cmd.forward()
        cmd.refresh()
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        cmd.meter_reset()      
        start = time.time()
        while elapsed<self.gl:
            cmd.frame(1)
            cmd.refresh()
            cmd.frame(2)
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('BLITS_V1',2*cnt/elapsed)

    def surface_calculation(self):
        cmd.load("$PYMOL_DATA/demo/il2.pdb")
        cmd.zoom(complete=1)
        cmd.hide()
        cmd.show("surface")
        cmd.clip("slab",0)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        start = time.time()
        while (elapsed)<self.short_cpu: 
            cmd.rebuild()
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('SURFACE_CALCULATION_V1',60*cnt/elapsed)

    def mesh_calculation(self):
        cmd.load("$PYMOL_DATA/demo/il2.pdb")
        cmd.zoom(complete=1)
        cmd.hide()
        cmd.show("mesh")
        cmd.clip("slab",0)
        cnt = 0
        elapsed = 0.0
        cmd.refresh()
        start = time.time()
        while (elapsed)<self.short_cpu: 
            cmd.rebuild()
            cmd.refresh()
            cnt = cnt + 1
            elapsed = time.time()-start
        report('MESH_CALCULATION_V1',60*cnt/elapsed)

    def ray_trace0(self): # Interactive benchmark
        self.configure()
        self.ray_tracing([
            [2,90],
            ])

    def ray_trace1(self): # Standard benchmark
        self.configure()
        self.ray_tracing([
            [1,90],
            [2,90],
            [4,90],
            [8,90],
            [1,120],
            [2,120],
            [1,160],
            [2,160],
            [1,200],
            [2,200],
            ])

    def ray_trace2(self): # Heavy-duty SMP workout
        self.configure()
        self.ray_tracing([
            [1,200],
            [2,200],
            [3,200],
            [4,200],
            [5,200],
            [6,200],
            [7,200],
            [8,200],
            [9,200],
            [10,200],
            [11,200],
            [12,200],
            ],width=3600,height=2700)
        
    def ray_tracing(self,conditions,width=640,height=480):
        cmd.load("$PYMOL_DATA/demo/1tii.pdb")
        cmd.zoom(complete=1)
        cmd.hide()
        cmd.show("spheres","11-15/")
        cmd.show("surface","21-25/")
        cmd.show("mesh","A/10-20/")
        cmd.show("sticks","41-50/")
        cmd.show("lines","51-55/")
        cmd.show("dots","61-65/")      
        cmd.show("cartoon","80-90/")
        cmd.turn('x',25)
        cmd.turn('y',25)
        for cond in conditions:
            (max_threads,hash_max) = cond
            cmd.set('max_threads',max_threads) 
            cmd.set('hash_max',hash_max) 
            cnt = 0
            elapsed = 0.0
            cmd.refresh()
            start = time.time()
            while elapsed<self.long_cpu:
                cmd.ray(width,height,quiet=1)
                cnt = cnt + 1
                elapsed = time.time()-start
            report('RAY_V2_PX%d_TH%02d_HSH%03d'%(width*height,
                                                                      max_threads,hash_max),60*cnt/elapsed)
        
    def get_prompt(self):
        self.prompt = self.message
        return self.prompt

    def delay_launch(self,action):
        self.configure()
        cmd.viewport(640,480)
        cmd.feedback("disable","all","everything")
        cmd.feedback("enable","python","output")
        t = threading.Thread(target=bench_fn,args=(action,))
        t.setDaemon(1)
        t.start()
        
    def get_panel(self):
        return [
            [ 1, 'Benchmarks', '' ],
            [ 2, 'Run All', 'cmd.get_wizard().delay_launch("run_all")' ],
            [ 2, 'Run GL', 'cmd.get_wizard().delay_launch("run_gl")' ],
            [ 2, 'Run CPU', 'cmd.get_wizard().delay_launch("run_cpu")' ],                  
            [ 2, 'Updates', 'cmd.get_wizard().delay_launch("updates")'],
            [ 2, 'Smooth Lines', 'cmd.get_wizard().delay_launch("smooth_lines")'],
            [ 2, 'Jagged Lines', 'cmd.get_wizard().delay_launch("jagged_lines")'],
            [ 2, 'Dots', 'cmd.get_wizard().delay_launch("dots")'],         
            [ 2, 'Jagged Lines', 'cmd.get_wizard().delay_launch("jagged_lines")'],
            [ 2, 'Sticks', 'cmd.get_wizard().delay_launch("sticks")'],
            [ 2, 'Surface', 'cmd.get_wizard().delay_launch("surface")'],
            [ 2, 'Spheres', 'cmd.get_wizard().delay_launch("spheres")'],
            [ 2, 'Cartoon', 'cmd.get_wizard().delay_launch("cartoon")'],
            [ 2, 'Blits', 'cmd.get_wizard().delay_launch("blits")'],                                    
            [ 2, 'Surface Calculation', 'cmd.get_wizard().delay_launch("surface_calculation")'],
            [ 2, 'Mesh Calculation', 'cmd.get_wizard().delay_launch("mesh_calculation")'],
            [ 2, 'Ray Tracing', 'cmd.get_wizard().delay_launch("ray_trace0")'],
            [ 2, 'End Demonstration', 'cmd.set_wizard()' ]
            ]

