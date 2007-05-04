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

import cmd
import math
import string
import pymol
import os
import traceback
from pymol import movie
# legacy mappings, remove in PyMOL 2.0

mload = movie.load
mrock = movie.rock
mroll = movie.roll

# should match the list in layer1/Color.c:
_color_cycle = [
    26   , # /* carbon */
    5    , # /* cyan */
    154  , # /* lightmagenta */
    6    , # /* yellow */
    9    , # /* salmon */
    29   , # /* hydrogen */
    11   , # /* slate */
    13   , # /* orange */
    10   , # /* lime */
    5262 , # /* deepteal */
    12   , # /* hotpink */
    36   , # /* yelloworange */
    5271 , # /* violetpurple */
    124  , # /* grey70 */
    17   , # /* marine */
    18   , # /* olive */
    5270 , # /* smudge */
    20   , # /* teal */
    5272 , # /* dirtyviolet */
    52   , # /* wheat */
    5258 , # /* deepsalmon */
    5274 , # /* lightpink */
    5257 , # /* aquamarine */
    5256 , # /* paleyellow */
    15   , # /* limegreen */
    5277 , # /* skyblue */
    5279 , # /* warmpink */
    5276 , # /* limon */
    53   , # /* violet */
    5278 , # /* bluewhite */
    5275 , # /* greencyan */
    5269 , # /* sand */
    22   , # /* forest */
    5266 , # /* lightteal */
    5280 , # /* darksalmon */
    5267 , # /* splitpea */
    5268 , # /* raspberry */
    104  , # /* grey50 */
    23   , # /* deepblue */
    51   , # /* brown */
    ]

_color_cycle_len = len(_color_cycle)

def mass_align(target,enabled_only=0,max_gap=50):
    list = cmd.get_names("public_objects",int(enabled_only))
    filter(lambda x:cmd.get_type(x)!="object:molecule",list)
    if enabled_only:
        aln_object = 'aln_enabled_to'+target
    else:
        aln_object = 'aln_all_to_'+target
    cmd.delete(aln_object)
    for name in list:
        if name!=target:
            if cmd.count_atoms("(%s) and (%s)"%(target,name))==0:
                cmd.align('polymer and name ca and (%s)'%name,
                'polymer and name ca and (%s)'%target,max_gap=max_gap,quiet=0,
                          object=aln_object)
    
def sum_formal_charges(selection="(all)",quiet=1):
    pymol.stored._util_sum_fc = 0.0
    cmd.iterate(selection,"stored._util_sum_fc=stored._util_sum_fc+formal_charge",quiet=1)
    result = pymol.stored._util_sum_fc
    if not quiet:
        print " util.sum_formal_charges: sum = %0.1f"%result
    return result

def sum_partial_charges(selection="(all)",quiet=1):
    pymol.stored._util_sum_pc = 0.0
    cmd.iterate(selection,"stored._util_sum_pc=stored._util_sum_pc+partial_charge",quiet=1)
    result = pymol.stored._util_sum_pc
    if not quiet:
        print " util.sum_partial_charges: sum = %0.4f"%result
    return result

def protein_assign_charges_and_radii(obj_name):

    from chempy.champ import assign

    # apply a few kludges
    
    # convent Seleno-methionine to methionine

    cmd.alter(obj_name+"///MSE/SE","elem='S';name='SD'",quiet=1)
    cmd.alter(obj_name+"///MSE/","resn='MET'",quiet=1)
    cmd.flag("ignore",obj_name,"clear")

    # remove alternate conformers

    cmd.remove(obj_name+" and not alt ''+A")
    cmd.alter(obj_name,"alt=''")
    cmd.sort(obj_name)
    cmd.fix_chemistry(obj_name,obj_name,1)
    
    # make sure all atoms are included...
    cmd.alter(obj_name,"q=1.0",quiet=1)
    
    print " Util: Fixing termini and assigning formal charges..."
    
    assign.missing_c_termini(obj_name,quiet=1)

    while not assign.formal_charges(obj_name,quiet=1):
        print " WARNING: unrecognized or incomplete residues are being deleted:"
        cmd.iterate("(byres ("+obj_name+" and flag 23)) and flag 31",
                        'print "  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/"',quiet=1)
        cmd.remove("byres ("+obj_name+" and flag 23)") # get rid of residues that weren't assigned
        assign.missing_c_termini(obj_name,quiet=1)
        
    print " Util: Assigning Amber 99 charges and radii..."
    
    cmd.h_add(obj_name)
    if not assign.amber99(obj_name,quiet=1):
        print " WARNING: some unassigned atoms are being deleted:"
        cmd.iterate("byres ("+obj_name+" and flag 23)",
                        'print "  "+model+"/"+segi+"/"+chain+"/"+resn+"`"+resi+"/"+name+"? ["+elem+"]"',quiet=1)
        cmd.remove(obj_name+" and flag 23") # get rid of any atoms that weren't assigned
        
    # show the user what the net charges are...
        
    formal = sum_formal_charges(obj_name,quiet=0)
    partial = sum_partial_charges(obj_name,quiet=0)
    if round(formal)!=round(partial):
        print " WARNING: formal and partial charge sums don't match -- there is a problem!"
    
def protein_vacuum_esp(selection, mode=2,border=10.0,quiet = 1):

    if ((string.split(selection)!=[selection]) or
         selection not in cmd.get_names('objects')):
        print " Error: must provide an object name"
        raise cmd.QuietException
    obj_name = selection + "_e_chg"
    map_name = selection + "_e_map"
    pot_name = selection + "_e_pot"
    cmd.disable(selection)
    cmd.delete(obj_name)
    cmd.delete(map_name)
    cmd.delete(pot_name)
    cmd.create(obj_name,"((polymer and ("+selection+
               ") and (not resn A+C+T+G+U)) or ((bymol (polymer and ("+
               selection+"))) and resn NME+NHE+ACE)) and (not hydro)")
         # try to just get protein...

    protein_assign_charges_and_radii(obj_name)
        
    ext = cmd.get_extent(obj_name)
    max_length = max(abs(ext[0][0] - ext[1][0]),abs(ext[0][1] - ext[1][1]),abs(ext[0][2]-ext[1][2])) + 2*border

    # compute an grid with a maximum dimension of 50, with 10 A borders around molecule, and a 1.0 A minimum grid

    sep = max_length/50.0
    if sep<1.0: sep = 1.0
    print " Util: Calculating electrostatic potential..."
    if mode==0: # absolute, no cutoff
        cmd.map_new(map_name,"coulomb",sep,obj_name,border)
    elif mode==1: # neutral, no cutoff
        cmd.map_new(map_name,"coulomb_neutral",sep,obj_name,border)
    else: # local, with cutoff
        cmd.map_new(map_name,"coulomb_local",sep,obj_name,border)      
        
    cmd.ramp_new(pot_name, map_name, selection=obj_name,zero=1)
    cmd.hide("everything",obj_name)
    cmd.show("surface",obj_name)
    cmd.set("surface_color",pot_name,obj_name)
    cmd.set("surface_ramp_above_mode",1,obj_name)
    
def color_carbon(color,selection="(all)"):
    selection = str(selection)
    cmd.color(color,"(%s) and elem c"%selection)
    
def cbss(selection="(all)",helix_color="red",sheet_color="yellow",loop_color="green",quiet=1):
    sel = str(selection)
    h = str(helix_color)
    s = str(sheet_color)
    l = str(loop_color)
    cmd.color(h,"(ss H and ("+sel+"))",quiet=quiet)
    cmd.color(s,"(ss S and ("+sel+"))",quiet=quiet)
    cmd.color(l,"((not (ss S+H)) and ("+sel+"))",quiet=quiet)

def cbag(selection="(all)",quiet=1):
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("carbon","(elem C and ("+s+"))",quiet=quiet)
    
def cbac(selection="(all)",quiet=1):
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("cyan","(elem C and ("+s+"))",quiet=quiet)
    
def cbam(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("lightmagenta","(elem C and ("+s+"))",quiet=quiet)

def cbay(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("yellow","(elem C and ("+s+"))",quiet=quiet)

def cbas(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("salmon","(elem C and ("+s+"))",quiet=quiet)

def cbaw(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("hydrogen","(elem C and ("+s+"))",quiet=quiet)

def cbab(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("slate","(elem C and ("+s+"))",quiet=quiet)

def cbao(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("brightorange","(elem C and ("+s+"))",quiet=quiet)

def cbap(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("purple","(elem C and ("+s+"))",quiet=quiet)

def cbak(selection="(all)",quiet=1):
    s = str(selection)   
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color("pink","(elem C and ("+s+"))",quiet=quiet)

def cnc(selection="(all)",quiet=1):
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)

def cba(color,selection="(all)",quiet=1):
    s = str(selection)
    cmd.color("atomic","(("+s+") and not elem C)",quiet=quiet)
    cmd.color(color,"(elem C and ("+s+"))",quiet=quiet)
    cmd.color(color,s,flags=1,quiet=quiet)

def performance(mode):
    mode = int(mode)
    if mode==0: # maximum quality
        cmd.set('line_smooth',1)
        cmd.set('depth_cue',1)         
        cmd.set('specular',1)
        cmd.set('surface_quality',1)
        cmd.set('stick_quality',15)
        cmd.set('sphere_quality',2)
        cmd.set('cartoon_sampling',14)
        cmd.set('ribbon_sampling',10)
        cmd.set('transparency_mode',2)
        cmd.set('stick_ball',1)
        cmd.do("rebuild")
    elif mode==33:
        cmd.set('line_smooth',1)         
        cmd.set('depth_cue',1)         
        cmd.set('specular',1)
        cmd.set('surface_quality',0)
        cmd.set('stick_quality',8)
        cmd.set('sphere_quality',1)
        cmd.set('cartoon_sampling',7)
        cmd.set('ribbon_sampling',1)
        cmd.set('transparency_mode',2)
        cmd.set('stick_ball',0)
        cmd.do("rebuild")
    elif mode==66: # good perfomance
        cmd.set('line_smooth',0)
        cmd.set('depth_cue',0)         
        cmd.set('specular',1)
        cmd.set('surface_quality',0)
        cmd.set('stick_quality',8)
        cmd.set('sphere_quality',1)
        cmd.set('cartoon_sampling',6)
        cmd.set('ribbon_sampling',1)
        cmd.set('transparency_mode',2)
        cmd.set('stick_ball',0.0)
        cmd.do("rebuild")         
    else: # maximum performance
        cmd.set('line_smooth',0)
        cmd.set('depth_cue',0)
        cmd.set('specular',0)
        cmd.set('surface_quality',-1) # new
        cmd.set('stick_quality',5)
        cmd.set('sphere_quality',0)
        cmd.set('ribbon_sampling',1)
        cmd.set('cartoon_sampling',3)
        cmd.set('transparency_mode',0)
        cmd.set('stick_ball',0.0)
        cmd.do("rebuild")         
    
    
def label_chains(sele="all"):
    pymol.stored._cs = []
    last = None
    save = ()
    list = []
    cmd.iterate(sele,"stored._cs.append((model,chain,index))")
    for a in pymol.stored._cs:
        if (a[0:2]!=save):
            list.append(last)
            list.append(a)
            save = a[0:2]
        last = a
    if len(list):
        list.append(last)
    list = filter(None,list)
    for a in list:
        if(a[1]==''):
            cmd.label("%s`%d"%(a[0],a[2]),'''"chain ''"''',quiet=1)         
        elif(a[1]==' '):
            cmd.label("%s`%d"%(a[0],a[2]),'''"chain ' '"''',quiet=1)         
        else:
            cmd.label("%s`%d"%(a[0],a[2]),"'chain '+chain",quiet=1)

def label_segments(sele="all"):
    pymol.stored._cs = []
    last = None
    save = ()
    list = []
    cmd.iterate(sele,"stored._cs.append((model,segi,index))")
    for a in pymol.stored._cs:
        if (a[0:2]!=save):
            list.append(last)
            list.append(a)
            save = a[0:2]
        last = a
    if len(list):
        list.append(last)
    list = filter(None,list)
    for a in list:
        if(a[1]==''):
            cmd.label("%s`%d"%(a[0],a[2]),'''"segi ''"''',quiet=1)         
        elif(a[1]==' '):
            cmd.label("%s`%d"%(a[0],a[2]),'''"segi ' '"''',quiet=1)         
        else:
            cmd.label("%s`%d"%(a[0],a[2]),"'segi '+segi",quiet=1)

    
def hide_sele():
    arg = cmd.get_names("selections")
    for a in arg:
        cmd.disable(a)

# FUBAR
#def hbond(a,b,cutoff=3.3,name='hbond'):
#   cmd.dist(name,"((%s) and ((%s) around %4.2f) and elem N,O)"%(a,b,cutoff),
#            "((%s) and ((%s) around %4.2f) and elem N,O)"%(b,a,cutoff),
#            cutoff)


def cbc(selection='(all)',first_color=7,quiet=1,legacy=0): 
    '''
    Color all chains a different color
    '''
    if int(legacy):
        c = first_color
        for a in cmd.get_chains(selection):
            if len(string.strip(a)):
                if not quiet: print (" util.cbc: color %d,(chain %s)"%(c,a))
                cmd.color("%d"%c,"(chain %s and (%s))"%(a,selection),quiet=quiet)
                c = c + 1
            elif len(a): # note, PyMOL's selection language can't handle this right now
                if not quiet: print (" util.cbc: color %d,(chain ' ')"%(c))
                cmd.color("%d"%c,"(chain '' and (%s))"%selection,quiet=quiet)
                c = c + 1
            else:
                if not quiet: print (" util.cbc: color %d,(chain '')"%(c))
                cmd.color("%d"%c,"(chain '' and (%s))"%selection,quiet=quiet)
                c = c + 1
    else:
        c = 0
        for a in cmd.get_chains(selection):
            if len(string.strip(a)):
                if not quiet: print (" util.cbc: color %d,(chain %s)"%(_color_cycle[c],a))
                cmd.color(_color_cycle[c],"(chain %s and (%s))"%(a,selection),quiet=quiet)
            elif len(a): # note, PyMOL's selection language can't handle this right now
                if not quiet: print (" util.cbc: color %d,(chain ' ')"%(_color_cycle[c]))
                cmd.color(_color_cycle[c],"(chain '' and (%s))"%selection,quiet=quiet)
            else:
                if not quiet: print (" util.cbc: color %d,(chain '')"%(_color_cycle[c]))
                cmd.color(_color_cycle[c],"(chain '' and (%s))"%selection,quiet=quiet)
            c = (c + 1) % _color_cycle_len
        
def color_objs(selection='(all)',quiet=1): 
    '''
    Color all chains a different color
    '''
    c = 0
    for a in cmd.get_names(selection=selection):
	if (selection!='all') and (selection!='(all)'):
            cmd.color(_color_cycle[c],"(?%s and (%s))"%(a,selection),quiet=quiet)
        else:
            cmd.color(_color_cycle[c],"(?%s)"%(a),quiet=quiet)
        c = (c + 1) % _color_cycle_len

def chainbow(selection='(all)',first_color=7): # NOT THREAD SAFE
    '''
    Color all chains in rainbow
    '''
    for a in cmd.get_chains(selection):
        if len(a):
            cmd.spectrum('count',selection="(chain %s and (%s))"%(a,selection),byres=1)
        else:
            cmd.spectrum('count',selection="(chain '' and (%s))"%selection,byres=1)
            
color_chains = cbc

def sum_charge(*arg,**kw): # NOT THREAD SAFE
    result = None
    try:
        obj = "all"
        if len(arg):
            obj = arg

        pymol.stored._sum_charge = 0.0
        cmd.iterate("(%s)"%obj,
                        "stored._sum_charge=stored._sum_charge+partial_charge")
        result = pymol.stored._sum_charge
        print " sum_charge: %6.4f"%result
    except:
        print " sum_charge: an error occurred."
    return result

def ray_shadows(mode):
    # adjustment factors for new lighting model in 0.99
    reflect_scale = 0.5
    direct_scale = 1.8
    gamma_scale = 1/1.3
    ambient_scale = 1.16
    if mode=='none':
        cmd.set('ray_shadows',0)
    else:
        cmd.set('ray_shadows',1)
    if mode=='light': 
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.65)
        cmd.set('reflect',0.25)
        cmd.set('shininess',40)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('spec_count',-1)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='matte':
        cmd.set('light_count',4)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set("light2","[-0.3,-0.4,-1.0]")
        cmd.set("light3","[-0.3,-0.3,-1.0]")
        cmd.set("light4","[-0.4,-0.3,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.45)
        cmd.set('reflect',0.45)
        cmd.set('shininess',25)
        cmd.set('power',1.25)
        cmd.set('spec_count',-1)                
        cmd.set('specular_intensity',0.2)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='soft':
        cmd.set('light_count',10)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set("light2","[-0.3,-0.4,-1.0]")
        cmd.set("light3","[-0.3,-0.3,-1.0]")
        cmd.set("light4","[-0.4,-0.3,-1.0]")
        cmd.set("light5","[-0.4,-0.5,-1.0]")
        cmd.set("light6","[-0.5,-0.5,-1.0]")
        cmd.set("light7","[-0.5,-0.4,-1.0]")
        cmd.set("light8","[-0.5,-0.3,-1.0]")
        cmd.set("light9","[-0.3,-0.5,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.40)
        cmd.set('reflect',0.50)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_count',1)
        cmd.set('shininess',55)
        cmd.set('power',1.0)
        cmd.set('spec_direct',0)        
        cmd.set('ray_shadow_decay_factor',0.1)
        cmd.set('ray_shadow_decay_range',1.8)
    elif mode=='occlusion':
        cmd.set('light_count',9)
        cmd.set("light" ,"[-0.2,-0.2,-1.0]")
        cmd.set("light2","[-0.2, 0.0,-1.0]")
        cmd.set("light3","[-0.2, 0.2,-1.0]")
        cmd.set("light4","[ 0.0, 0.2,-1.0]")
        cmd.set("light5","[ 0.2, 0.2,-1.0]")
        cmd.set("light6","[ 0.2, 0.0,-1.0]")
        cmd.set("light7","[ 0.2,-0.2,-1.0]")
        cmd.set("light8","[ 0.0,-0.2,-1.0]")
        cmd.set('ambient',0.18)
        cmd.set('direct',0.10)
        cmd.set('reflect',0.80)
        cmd.set('shininess',10)
        cmd.set('spec_count',-1)        
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0)
        cmd.set('spec_direct',0.25)
        cmd.set('ray_shadow_decay_factor',0.1)
        cmd.set('ray_shadow_decay_range',5.0)
    elif mode=='medium':
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.14)
        cmd.set('direct',0.45)
        cmd.set('reflect',0.45)
        cmd.set('shininess',55)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='heavy':
        cmd.set('light_count',2)        
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.05)
        cmd.set('direct',0.20)
        cmd.set('reflect',0.85)
        cmd.set('spec_count',-1)
        cmd.set('shininess',90) 
        cmd.set('power',1.0) 
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
    elif mode=='black': # best for light backgrounds
        cmd.set('light_count',2)
        cmd.set("light","[-0.4,-0.4,-1.0]")
        cmd.set('ambient',0.001)
        cmd.set('direct',0.0)
        cmd.set('reflect',1.1)
        cmd.set('spec_count',-1)
        cmd.set('power',1.0)
        cmd.set('shininess',90)
        cmd.set('specular_intensity',0.5)
        cmd.set('spec_direct',0)
        cmd.set('ray_shadow_decay_factor',0)
        
def ff_copy(src,dst): # NOT THREAD SAFE
    pymol._rcopy = pymol.Scratch_Storage()
    pymol._rcopy.pc={}
    pymol._rcopy.tt={}
    cmd.iterate("(%s)"%src,"_rcopy.pc[name]=partial_charge")
    cmd.alter("(%s)"%dst,"partial_charge=_rcopy.pc[name]")
    cmd.iterate("(%s)"%src,"_rcopy.tt[name]=text_type")
    cmd.alter("(%s)"%dst,"text_type=_rcopy.tt[name]")
    del pymol._rcopy
    
def b2vdw(*arg,**kw):
    if not len(arg):
        sele = 'all'
    else:
        sele = arg[0]
    # use B values to create RMS VDW spheres
    # rms = sqrt(b/(8*(PI^2)))
    cmd.alter("(%s)"%sele,"vdw=math.sqrt(b/78.9568352087)")
    
def phipsi(selection="(pk1)"): # NOT THREAD SAFE
    n_sele =   "((byres (%s)) & name n)"%selection
    c_sele =   "((byres (%s)) & name c)"%selection
    ca_sele =  "((byres (%s)) & name ca)"%selection
    cm_sele = "((neighbor (%s)) and not (byres (%s)))"%(n_sele,n_sele)
    np_sele = "((neighbor (%s)) and not (byres (%s)))"%(c_sele,c_sele)
    cmd.feedback("push")
    cmd.feedback("disable","selector","everythin")
    cm_cnt = cmd.select("_pp_cm",cm_sele)
    n_cnt = cmd.select("_pp_n",n_sele)
    c_cnt = cmd.select("_pp_c",c_sele)
    ca_cnt = cmd.select("_pp_ca",ca_sele)
    np_cnt = cmd.select("_pp_np",np_sele)
    if(cm_cnt and n_cnt and ca_cnt and c_cnt):
        phi = cmd.get_dihedral("_pp_c","_pp_ca","_pp_n","_pp_cm")
    else:
        phi = None
    if(n_cnt and ca_cnt and c_cnt and np_cnt):
        psi = cmd.get_dihedral("_pp_np","_pp_c","_pp_ca","_pp_n")
    else:
        psi = None
    cmd.feedback("pop")
    cmd.delete("_pp_cm")
    cmd.delete("_pp_n")
    cmd.delete("_pp_c")
    cmd.delete("_pp_ca")
    cmd.delete("_pp_np")
    return (phi,psi)

def rainbow(selection="(name ca and alt '',A)",reverse=0): # NOT THREAD SAFE

    cmd.feedback("push")
    cmd.feedback("disable","executive","actions")

    # your basic rainbow...
    
    list = [
        (0,0,255),      
        (0,0,255),
        (0,128,255),
        (0,255,255),
        (0,255,128),           
        (0,255,0),
        (128,255,0),
        (255,255,0),
        (255,128,0),
        (255,0,0),
        (255,0,0)      
        ]
    if reverse:
        list.reverse()
    #
    last = list.pop(0)
    cmd.set_color("_000",[last[0]/255.0,last[1]/255.0,last[2]/255.0])
    c = 1
    for a in list:
        for b in range(1,21):
            b0 = b/20.0
            b1 = 1.0-b0
            cname = "_%03d"%c
            r = last[0]*b1+a[0]*b0
            g = last[1]*b1+a[1]*b0
            b = last[2]*b1+a[2]*b0
            cmd.set_color(cname,[r/255.0,g/255.0,b/255.0])
            c = c + 1
        last = a

    cas = cmd.index("((byres ("+selection+")) and name ca and not het)")
    l = len(cas)
    if not len(cas):
        return
    c = 0
    for a in cas:
        col = int((200*c)/l)
        cmd.color("_%03d"%col,"((%s) and (byres %s`%d))"%(selection,a[0],a[1]))
        c = c + 1

    cmd.feedback("pop")
    
def ss(selection="(name ca and alt '',A)",state=1): # NOT THREAD SAFE

    print ' util.ss: WARNING: This is not a "correct" secondary structure'
    print ' util.ss: assignment algorithm!  Please use only as a last resort.'
    
    cmd.feedback("push")
    cmd.feedback("disable","executive","actions")
    
    ss_pref = "_sss"
    sss1 = ss_pref+"1"
    cnt = cmd.select(sss1,"((byres ("+selection+")) and name ca and not het)")
    print " util.ss: initiating secondary structure assignment on %d residues."%cnt
    cas = cmd.index(sss1)
    if not len(cas):
        return
    # set cartoon mode to auto over the selection
    
    cmd.cartoon("auto",sss1)

    print " util.ss: extracting sequence and relationships..."

    # get CA list
    
    res_list = []
    pymol._ss = pymol.Scratch_Storage()
    pymol._ss.res_list = res_list
    cmd.iterate(sss1,'_ss.res_list.append((model,index))')

    # generate atom-to-residue conversion dictionaries

    ca_dict = {}
    n_dict = {}
    o_dict = {}
    scr_dict = {} # scr = segment,chain,resi 
    pymol._ss.n_dict = n_dict
    pymol._ss.o_dict = o_dict
    pymol._ss.scr_dict = scr_dict
    pymol._ss.ca_dict = ca_dict
    cmd.iterate(sss1,
                    '_ss.scr_dict[(model,index)]=(segi,chain,resi)') # CA's
    cmd.iterate("((byres "+sss1+") and n;n)"
                    ,'_ss.scr_dict[(model,index)]=(segi,chain,resi)') # N's
    cmd.iterate("((byres "+sss1+") and n;o)",
                    '_ss.scr_dict[(model,index)]=(segi,chain,resi)') # O's
    cmd.iterate(sss1,
                    '_ss.ca_dict[(segi,chain,resi)] = (model,index)')
    cmd.iterate("((byres "+sss1+") and n;n)",
                    '_ss.n_dict[(segi,chain,resi)] = (model,index)')
    cmd.iterate("((byres "+sss1+") and n;o)",
                    '_ss.o_dict[(segi,chain,resi)] = (model,index)')

    scr_dict[None]=None
    o_dict[None]=None
    n_dict[None]=None
    ca_dict[None]=None
    
    # create special version of cas with gaps

    gap = [None,None,None,None]  
    # gap large enough to distinguish i+4 interations from gaps
    last = None
    for a in res_list:
        if last!=None:
            if(cmd.count_atoms(
                "((neighbor(neighbor(neighbor (%s`%d)))) and (%s`%d))"%
                (last[0],last[1],a[0],a[1]),quiet=1)==0):
                gap.extend([None,None,None,None])
        gap.append(a)
        last = a
    gap.extend([None,None,None,None])

    print " util.ss: analyzing phi/psi angles (slow)..."

    # generate reverse-lookup for gap indices

    ss = {}

    c = 0
    gap_idx = {}
    for a in gap:
        gap_idx[a] = c
        c = c + 1

    # secondary structure database...
    
    ss = {}
    ss[None]=None
    
    # make decisions based on phi/psi

    for a in cas:
        ss[a] = 'L' # default
    phipsi = cmd.get_phipsi(sss1,state)
    for a in phipsi.keys():
        (phi,psi) = phipsi[a]
#      print scr_dict[a],(phi,psi)
        if (phi!=None) and (psi!=None):
            if ((phi<-45) and (phi>-160) and
                 (psi<-170) or (psi>10)): # beta?
                ss[a] = 's'
            elif ((phi<-45) and (phi>-160) and
                    (psi>-80) and (psi<-25)): # helix?
                ss[a] = 'H'
                
    print " util.ss: finding hydrogen bonds..."
    
    # find all pairwise hydrogen bonds and make note of them in dict

    hb = cmd.find_pairs("((byres "+sss1+") and n;n)",
                              "((byres "+sss1+") and n;o)",mode=1,
                              cutoff=3.7,angle=55,
                              state1=state,state2=state)
    
    hb_dict = {}  # [((N-atom) (O-atom))] = 1
    n_hb_dict = {} # [(N-atom)] = [(O-atom),...]
    o_hb_dict = {} # [(O-atom)] = [(N-atom),...]
    for a in hb:
#      cmd.dist("(%s`%d)"%a[0],"(%s`%d)"%a[1])
        hb_dict[a] = 1
        n = a[0]
        o = a[1]
        if not n_hb_dict.has_key(n): n_hb_dict[n]=[]
        if not o_hb_dict.has_key(o): o_hb_dict[o]=[]
        n_hb_dict[n].append(o)
        o_hb_dict[o].append(n)

    # check to insure that all helical residues have at least an i +/- 4
    # hydrogen bond

    for c in xrange(4,len(gap)-4):
        a = gap[c]
        if ss[a]=='H':
            aN = n_dict[scr_dict[a]]
            aO = o_dict[scr_dict[a]]
            am4O = o_dict[scr_dict[gap[c-4]]]
            ap4N = n_dict[scr_dict[gap[c+4]]]
            if not hb_dict.has_key((aN,am4O)):
                if not hb_dict.has_key((ap4N,aO)):
                    ss[a]='L'

    print " util.ss: verifying beta sheets..."
    
    # check to insure that all beta residues have proper interactions

    rep_dict = {}
    repeat = 1
    while repeat:
        repeat = 0
        c = 4
        cc = len(gap)-4
        while c<cc:
            a1 = gap[c]
            if (ss[a1] in ['s','S']) and not rep_dict.has_key(a1):
                rep_dict[a1] = 1
                valid = 0
                scr_a1 = scr_dict[a1]
                # look for antiparallel 2:2 H-bonds (NH-O=C + C=O-HN) 
                n_a1_atom = n_dict[scr_a1]
                o_a1_atom = o_dict[scr_a1]
                if (n_hb_dict.has_key(n_a1_atom) and 
                     o_hb_dict.has_key(o_a1_atom)):
                    for n_hb_atom in n_hb_dict[n_a1_atom]:
                        for o_hb_atom in o_hb_dict[o_a1_atom]:
                            n_hb_scr = scr_dict[n_hb_atom]
                            o_hb_scr = scr_dict[o_hb_atom]
                            if o_hb_scr == n_hb_scr:
                                b1 = ca_dict[o_hb_scr]
                                if abs(c-gap_idx[b1])>2:
                                    ss[b1] = 'S' 
                                    ss[a1] = 'S' 
                                    valid = 1
                # look for antiparallel offset HB (i,i+2,j,j-2)
                a3 = gap[c+2]
                if (a3!=None):
                    scr_a3 = scr_dict[a3]
                    o_a1_atom = o_dict[scr_a1]
                    n_a3_atom = n_dict[scr_a3]
                    if (n_hb_dict.has_key(n_a3_atom) and
                         o_hb_dict.has_key(o_a1_atom)):               
                        for n_hb_atom in n_hb_dict[n_a3_atom]:
                            for o_hb_atom in o_hb_dict[o_a1_atom]:
                                n_hb_scr = scr_dict[n_hb_atom]
                                o_hb_scr = scr_dict[o_hb_atom]
                                b1 = ca_dict[o_hb_scr]
                                if b1!=None:
                                    b1_i = gap_idx[b1]
                                    if abs(c-b1_i)>2: # no turns!
                                        b3 = gap[b1_i-2]
                                        if b3!=None:
                                            b3_scr = scr_dict[b3]
                                            if b3_scr == n_hb_scr:
                                                a2 = gap[c+1]
                                                b2 = gap[gap_idx[b1]-1]
                                                ss[b1] = 'S'
                                                ss[b3] = 'S'
                                                ss[a1] = 'S'
                                                ss[a3] = 'S'
                                                if ss[a2]=='L': ss[a2] = 's'
                                                if ss[b2]=='L': ss[b2] = 's'
                                                valid = 1
                # look for antiparallel offset HB (i,i-2,j,j+2)
                a3 = gap[c-2]
                if (a3!=None):
                    scr_a3 = scr_dict[a3]
                    n_a1_atom = n_dict[scr_a1]
                    o_a3_atom = o_dict[scr_a3]
                    if (n_hb_dict.has_key(n_a1_atom) and
                         o_hb_dict.has_key(o_a3_atom)):               
                        for n_hb_atom in n_hb_dict[n_a1_atom]:
                            for o_hb_atom in o_hb_dict[o_a3_atom]:
                                n_hb_scr = scr_dict[n_hb_atom]
                                o_hb_scr = scr_dict[o_hb_atom]
                                b1 = ca_dict[o_hb_scr]
                                if b1!=None:
                                    b1_i = gap_idx[b1]
                                    if abs(c-b1_i)>2: # no turns!
                                        b3 = gap[b1_i-2]
                                        if b3!=None:
                                            b3_scr = scr_dict[b3]
                                            if b3_scr == n_hb_scr:
                                                a2 = gap[c-1]
                                                b2 = gap[gap_idx[b1]-1]
                                                ss[b1] = 'S'
                                                ss[b3] = 'S'
                                                ss[a1] = 'S'
                                                ss[a3] = 'S'
                                                if ss[a2]=='L': ss[a2] = 's'
                                                if ss[b2]=='L': ss[b2] = 's'
                                                valid = 1
                # look for parallel 1:3 HB (i,j-1,j+1)
                n_a1_atom = n_dict[scr_a1]
                o_a1_atom = o_dict[scr_a1]
                if (n_hb_dict.has_key(n_a1_atom) and
                     o_hb_dict.has_key(o_a1_atom)):
                    for n_hb_atom in n_hb_dict[n_a1_atom]:
                        for o_hb_atom in o_hb_dict[o_a1_atom]:
                            n_hb_scr = scr_dict[n_hb_atom]
                            o_hb_scr = scr_dict[o_hb_atom]
                            b0 = ca_dict[n_hb_scr]
                            if b0!=None:
                                b2 = gap[gap_idx[b0]+2]
                                if b2!=None:
                                    b2_scr = scr_dict[b2]
                                    if b2_scr == o_hb_scr:
                                        b1 = gap[gap_idx[b0]+1]
                                        ss[a1] = 'S' 
                                        ss[b0] = 'S'
                                        if ss[b1]=='L': ss[b1]='s'
                                        ss[b2] = 'S'
                                        valid = 1
                                        repeat = 1
                if not valid:
                    ss[a1] = 'L'
            c = c + 1

    # automatically fill 1 residue gaps in helices and well-defined sheets
    c = 4
    cc = len(gap)-6
    while c<cc:
        a1 = gap[c]
        a3 = gap[c+2]
        ss_a1 = ss[a1]
        ss_a3 = ss[a3]
        if (ss_a1==ss_a3) and (ss_a1 in ['S','H']):
            a2 = gap[c+1]
            ss[a2] = ss_a1
        c = c + 1

    # remove singleton sheet residues
    c = 4
    cc = len(gap)-4
    while c<cc:
        a0 = gap[c-1]
        a1 = gap[c]
        a2 = gap[c+1]
        if ss[a1] in ['s','S']:
            if ((not ss[a0] in ['s','S']) and
                 (not ss[a2] in ['s','S'])):
                 ss[a1] = 'L'
        c = c + 1

    # remove sheet residues which aren't next to another sheet 
    c = 4
    cc = len(gap)-4
    while c<cc:
        a1 = gap[c]
        if ss[a1]=='S':
            a1 = gap[c]
            scr_a1 = scr_dict[a1]
            # look for hydrogen bonds to another sheet
            n_a1_atom = n_dict[scr_a1]
            o_a1_atom = o_dict[scr_a1]
            certain = 0
            if n_hb_dict.has_key(n_a1_atom):
                for n_hb_atom in n_hb_dict[n_a1_atom]:
                    n_hb_ca_atom=ca_dict[scr_dict[n_hb_atom]]
                    if ss[n_hb_ca_atom]=='S':
                        certain = 1
                        break
            if o_hb_dict.has_key(o_a1_atom):
                for o_hb_atom in o_hb_dict[o_a1_atom]:
                    o_hb_ca_atom=ca_dict[scr_dict[o_hb_atom]]
                    if ss[o_hb_ca_atom]=='S':
                        certain = 1
                        break
            if not certain:
                ss[a1] = 's'
        c = c + 1

    # remove questionable sheet residues
    c = 4
    cc = len(gap)-4
    while c<cc:
        a0 = gap[c-1]
        a1 = gap[c]
        a2 = gap[c+1]
        if ss[a1]=='s':
            if (not ((ss[a0]=='S') and (ss[a2]=='S'))):
                ss[a1] = 'L'
        c = c + 1

    # extend helices if hydrogen bonding requirements are met
    rep_dict = {}
    repeat = 1
    while repeat:
        repeat = 0
        c = 4
        cc = len(gap)-4
        while c<cc:
            a = gap[c]
            if not rep_dict.has_key(a):
                if ss[gap[c+1]]=='H':
                    rep_dict[a] = 1
                    if ss[a]!='H': # N-terminal end
                        aO = o_dict[scr_dict[a]]
                        ap4N = n_dict[scr_dict[gap[c+4]]]
                        ap3N = n_dict[scr_dict[gap[c+3]]]
                        if hb_dict.has_key((ap4N,aO)) or hb_dict.has_key((ap3N,aO)):
                            ss[a]='H'
                            repeat = 1
                            c = c - 5
                            if c<4: c=4
                if ss[gap[c-1]]=='H':
                    a = gap[c]
                    if ss[a]!='H': # C-terminal end
                        rep_dict[a] = 1
                        aN = n_dict[scr_dict[a]]
                        am4O = o_dict[scr_dict[gap[c-4]]]
                        am3O = o_dict[scr_dict[gap[c-3]]]
                        if hb_dict.has_key((aN,am4O)) or hb_dict.has_key((aN,am3O)):
                            ss[a]='H'
                            repeat = 1
                            c = c - 5
                            if c<4: c=4
            c = c + 1

    # remove doubleton helices

    c = 4
    cc = len(gap)-5
    while c<cc:
        a0 = gap[c-1]
        a1 = gap[c]
        a2 = gap[c+1]
        a3 = gap[c+2]
        ss_a0 = ss[gap[c-1]]
        ss_a1 = ss[gap[c]]
        ss_a2 = ss[gap[c+1]]
        ss_a3 = ss[gap[c+2]]
        if ss_a1=='H':
            if (ss_a2==ss_a1) and (ss_a0!=ss_a2) and (ss_a2!=ss_a3):
                ss[a1] = 'L'
                ss[a2] = 'L'
        c = c + 1

    # remove totally unreasonable helix and sheet residues

    c = 4
    cc = len(gap)-5
    while c<cc:
        a1 = gap[c]
        ss_a1 = ss[gap[c]]
        if ss_a1=='H':
            if phipsi.has_key(a1):
                (phi,psi) = phipsi[a1]
                if (phi>0) and (phi<150):
                    ss[a1] = 'L'
                elif((psi<-120) or (psi>140)):
                    ss[a1] = 'L'
        elif ss_a1 in ['S','s']:
            if phipsi.has_key(a1):
                (phi,psi) = phipsi[a1]
                if (phi>45) and (phi<160):
                    ss[a1] = 'L'
#            if (psi<-30) and (psi>-150):
                if (psi<-65) and (psi>-150):
                    ss[a1] = 'L'
            
        c = c + 1


    for x in range(1,3):
        # remove singleton sheet residues
        c = 4
        cc = len(gap)-4
        while c<cc:
            a0 = gap[c-1]
            a1 = gap[c]
            a2 = gap[c+1]
            if ss[a1] in ['s','S']:
                if ((not ss[a0] in ['s','S']) and
                     (not ss[a2] in ['s','S'])):
                     ss[a1] = 'L'
            c = c + 1

        # remove sheet residues which aren't next to another sheet 
        c = 4
        cc = len(gap)-4
        while c<cc:
            a1 = gap[c]
            if ss[a1]=='S':
                a1 = gap[c]
                scr_a1 = scr_dict[a1]
                # look for hydrogen bonds to another sheet
                n_a1_atom = n_dict[scr_a1]
                o_a1_atom = o_dict[scr_a1]
                certain = 0
                if n_hb_dict.has_key(n_a1_atom):
                    for n_hb_atom in n_hb_dict[n_a1_atom]:
                        n_hb_ca_atom=ca_dict[scr_dict[n_hb_atom]]
                        if ss[n_hb_ca_atom]=='S':
                            certain = 1
                            break
                if o_hb_dict.has_key(o_a1_atom):
                    for o_hb_atom in o_hb_dict[o_a1_atom]:
                        o_hb_ca_atom=ca_dict[scr_dict[o_hb_atom]]
                        if ss[o_hb_ca_atom]=='S':
                            certain = 1
                            break
                if not certain:
                    ss[a1] = 's'
            c = c + 1

        # remove questionable sheet residues
        c = 4
        cc = len(gap)-4
        while c<cc:
            a0 = gap[c-1]
            a1 = gap[c]
            a2 = gap[c+1]
            if ss[a1]=='s':
                if (not ((ss[a0]=='S') and (ss[a2]=='S'))):
                    ss[a1] = 'L'
            c = c + 1

#      lst = ss.keys()
#      lst.sort()
#      for a in lst: print scr_dict[a],ss[a]
        
    # assign protein
    for a in cas:
        if ss[a]=='s':
            ss[a]='S'
        
    cmd.alter(sss1,"ss ='L'")
    for a in cas:
        if ss[a]!='L':
            cmd.alter("(%s`%d)"%a,"ss='%s'"%ss[a])

    cmd.feedback("pop")

    del pymol._ss # IMPORTANT
    cmd.delete(sss1)
    cmd.rebuild(selection,'cartoon')
    #
#   print conn_hash.keys()
    print " util.ss: assignment complete."

def colors(scheme=""):
    if scheme=="jmol":
        cmd.set("auto_color",0)
        cmd.set_color("hydrogen",[1.000,1.000,1.000])
        cmd.set_color("carbon",[0.567,0.567,0.567])
        cmd.set_color("nitrogen",[0.189,0.315,0.976])
        cmd.set_color("oxygen",[1.000,0.051,0.051])
        cmd.set_color("fluorine",[0.567,0.882,0.314])
        cmd.set_color("sulfur",[1.000,1.000,0.189])
        cmd.color("carbon","elem c")
        cmd.recolor()
        
                          
