#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2003 by Warren Lyford Delano of DeLano Scientific. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* Scott Dixon, Metaphorics LLC
#-* 
#-*
#Z* -------------------------------------------------------------------

# For code contributed by Scott Dixon, the following notice applies:
# This source code is contributed to the public domain and may be freely
# copied and distributed for research, profit, fun or any other reason,
# with these restrictions: (1) unmodified or functionally equivalent code
# derived from this code must contain this notice, (2) all derived code
# must acknowledge the author and institution, and (3) no liability is
# assumed by the author(s) for any use or misuse of this software.

import cmd
from cmd import _cmd, is_tuple

from chempy import Storage,Atom,Bond
from chempy.models import Indexed
from pymol import util

import re
import string

from chempy import cex
CEX=cex

class CEXpyParser(CEX.CEXsmilesParser): # Author: Scott Dixon
    def __init__(self):
        self.model = Indexed()
        self.atomN = 0
    def MakeAtom(self, atnum):
        at = Atom()
        at.index = self.atomN
        self.atomN = self.atomN+1
        at.symbol = CEX.CEXsmilesParser.num2sym(self,atnum)
        self.model.atom.append(at)
        return at
    def MakeBond(self,at1, at2, bo):
        bnd = Bond()
        bnd.index = [ at1.index, at2.index ]
        bnd.order = bo
        self.model.bond.append(bnd)
    def SetFormalCharge(self, atom, charge):
        atom.formal_charge = charge
    def SetHcount(self, atom, count):
        pass
    def SetAtomicMass(self, atom, mass):
        pass

   
#---------------------------------------------------------------------------------
def readcex(file,*args):  # Author: Scott Dixon
    import os.path

      
    standardResidues = ("ALA", "ARG", "ASP", "ASN", "ASX", "CYS", "GLY", "GLU",
                        "GLN", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
                        "SER", "THR", "TYR", "TRP", "VAL", "HID", "HIE")
    def childprop(tag, tree,_self=cmd):
        try:
            p = CEX.selectChildren(tree,tag)[0].value
            if p[0] == '"': p = p[1:-1]
            list = map(string.strip,string.split(p,";"))
            return list
        except IndexError:
            return None
        
    f = open(file,"r")
    cs = CEX.CEXstream(f)
    modelname = os.path.splitext(os.path.basename(file))[0] #get the base filename
    if len(args) > 0: modelname = args[0]
    while 1:
        tree = CEX.readTree(cs)
        if not tree: break
        if tree.name == "$MOL":
            smiles = tree.value
            confs = CEX.selectChildren(tree,"XYZ")
            try:
                name = CEX.selectProperty(tree,"/NAM").value
            except AttributeError:
                name = ""
            anames = childprop("ANAM", tree)
            atypes = childprop("ATYPE", tree)
            resnos = childprop("RESNO", tree)
            resnames = childprop("RESNA", tree)
            chains = childprop("CHAIN", tree)
            occups = childprop("OCCUP", tree)
            bvals = childprop("BVAL", tree)
            pcharges = childprop("PCHARGE", tree)
            rads = childprop("RAD", tree)
            ligsites = childprop("SCORE", tree)
            clusters = childprop("CLUSTER", tree)
            parser = CEXpyParser()
            try:
                parser.parse(smiles)
            except CEX.CEXsmilesError, data:
                print data
                break
            parser.model.molecule.title = name
            for at in parser.model.atom:
                i = at.index
                if anames: at.name = anames[i]
                if resnames:
                    at.resn = resnames[i]
                    if at.resn in standardResidues: at.hetatm = 0
                if resnos:
                    at.resi = resnos[i]
                    at.resi_number = int(resnos[i])
                if chains: at.chain = chains[i]
                if bvals: at.b = bvals[i]
                if occups: at.q = occups[i]
                if rads: at.vdw = rads[i]
                if atypes: at.text_type = atypes[i]
                if pcharges:
                    chg = pcharges[i]
                    if chg == "": chg = 0.0
                    at.partial_charge = float(chg)
                if clusters:
                    at.resi = clusters[i]
                    at.resi_number = int(clusters[i])
                if ligsites: at.b = ligsites[i]
            for conf in confs:
                if conf.value[0] == '"':
                    conf.value = conf.value[1:-1]
                    xyz = map(lambda x:string.split(x,","),string.split(conf.value,";"))
                    for at in parser.model.atom:
                        i = at.index
                        at.coord = [ float(xyz[i][0]), float(xyz[i][1]), float(xyz[i][2]) ]
                cmd.load_model(parser.model,modelname,discrete=1,finish=0)
        else:
            pass
    _cmd.finish_object(_self._COb,modelname)
    f.close()
 
def colorbyB(selection="spheres",first=7,last=3): # Author: Scott Dixon
    cols = [(0.,0.,1.),(.5,.5,1.),(.8,.8,1.),(1.,1.,1.),
            (1.,.9,.9),(1.,.6,.6),(1.,0.,0.)]
    nbins = len(cols)
    incr = float(first-last)/nbins
    for i in range(nbins):
        cname = "_spcol%03d"%i
        cmd.set_color(cname,cols[i])
        b = first - i*incr
        cmd.color(cname,"((%s) and not(b > %f) and b > %f)"%(selection,b,b-incr))

def metaphorics(): 
    cmd.extend("readcex",readcex)
    cmd.extend("colorbyB",colorbyB)

def get_context_info():  # Author: Warren DeLano
    context_dict= {}
    context_list= []
    for a in cmd.get_names("all"):
       context = None   
       if a[-6:]=='_water': 
          context = a[:-6] 
       if a[-7:]=='_ligand':  
          context = a[:-7] 
       if a[-5:]=='_site':  
          context = a[:-5] 
       if a[-6:]=='_hbond':  
          context = a[:-6]
       if context!=None:
           if not context_dict.has_key(context):
               context_list.append(context)
               context_dict[context] = [a]
           else:
               context_dict[context].append(a)
    return (context_list,context_dict)

def toggle_labels(mode=-1):
    global labels
    if mode<0:
        labels = not labels
    else:
        labels = mode
    if labels:
        cmd.show("labels")
    else:
        cmd.hide("labels")

def toggle_cgos(mode=-1):
    global cgos
    if mode<0:
        cgos = not cgos
    else:
        cgos = mode
    if cgos:
        cmd.show("cgo")
    else:
        cmd.hide("cgo")

def toggle_ligands(mode=-1):
    global ligands
    if mode<0:
        ligands = not ligands
    else:
        ligands = mode
    if ligands:
        cmd.show("sticks",m4x_ligands)
    else:
        cmd.hide("sticks",m4x_ligands)

def toggle_sites(mode=-1):
    global sites
    if mode<0:
        sites = not sites
    else:
        sites = mode
    if sites:
        cmd.show("lines",m4x_sites)
    else:
        cmd.hide("lines",m4x_sites)

def toggle_waters(mode=-1):
    global waters
    if mode<0:
        waters = not waters
    else:
        waters = mode
    if waters:
        cmd.show("nonbonded",m4x_waters)
    else:
        cmd.hide("nonbonded",m4x_waters)

def toggle_dashes(mode=-1):
    global dashes
    if mode<0:
        dashes = not dashes
    else:
        dashes = mode
    if dashes:
        cmd.show("dashes")
    else:
        cmd.hide("dashes")

def toggle_zooms(mode=-1):
    global zooms
    if mode<0:
        zooms = not zooms
    else:
        zooms = mode
    if zooms:
        cmd.zoom(m4x_ligands,2)
    else:
        cmd.zoom("m4x_aligned")

def setup_contexts(context_info):   # Author: Warren DeLano
    (list,dict) = context_info[0:2]
    key_list = [
        'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10', #,'F11','F12',
        'SHFT-F1','SHFT-F2','SHFT-F3','SHFT-F4','SHFT-F5','SHFT-F6','SHFT-F7',
        'SHFT-F8','SHFT-F9','SHFT-F10']# ,'SHFT-F11','SHFT-F12']
    doc_list = ["Keys"]
    zoom_context = 1
    global labels
    labels = 1
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_labels)
        doc_list.append(key+": Toggle Dist")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,lambda :(cmd.zoom(),toggle_labels(0)))
        doc_list.append(key+": Zoom All")
    
    for a in list:
        water = a+"_water"
        ligand = a+"_ligand"
        site = a+"_site"
        hbond = a+"_hbond"
        name_list = dict[a]
        zoom_list = []
        if water in name_list:
            cmd.show("nonbonded",water)
            util.cbac(water)
            zoom_list.append(water)
        if ligand in name_list:
            cmd.show("sticks",ligand)
            cmd.hide("cartoon",ligand)
            util.cbag(ligand)
            zoom_list.append(ligand)
        if site in name_list:
            cmd.show("sticks",site)
            util.cbac(site)
            zoom_list.append(site)
            # replace cartoon with explicit atoms for "site" atoms
            cmd.hide("cartoon",site)
            cmd.show("sticks","(byres (neighbor ("+site+" and name c))) and name n+ca")
            cmd.show("sticks","(byres (neighbor ("+site+" and name n))) and name c+ca+o")
        if len(zoom_list):
            if len(key_list):
                key = key_list.pop(0)
                zoom_str = string.join(zoom_list,' or ')
                if zoom_context == 1:
                    zoom_context = zoom_str
                elif zoom_context not in (0,1):
                    zoom_context = 0
                cmd.set_key(key,lambda x=zoom_str:(cmd.zoom(x)))
                mo = re.search("_([^_]+)$",a)
                if mo:
                    cont_name = mo.groups()[0]
                else:
                    cont_name = a
                doc_list.append(key+": Zoom "+cont_name)
                
        if hbond in name_list:
            cmd.show("dashes",hbond)
            cmd.show("labels",hbond)

        
    cmd.wizard("fedora",doc_list)
    if zoom_context not in (0,1):
        cmd.zoom(zoom_context)
    toggle_labels(0)
#    cmd.feedback("enable","python","output")
    cmd.feedback("enable","objectmolecule","results")
    cmd.feedback("disable","selector","actions")
    cmd.feedback("disable","scene","actions")
    cmd.set("internal_feedback",1)
    cmd.set("internal_prompt",0)
    
        
def setup_alignment_contexts(context_info):   # Author: Warren DeLano
    (list,dict) = context_info[0:2]
    doc_list = ['\888Legend:']
    obj_name_dict = {}
    for a in list:
        sf = string.find(a,"_")
        if sf>=0:
            object_name = a[0:sf]
            if not obj_name_dict.has_key(object_name):
                obj_name_dict[object_name] = 1
                col_index = cmd.get_object_color_index(object_name)
                if col_index>=0:
                    col_tup = cmd.get_color_tuple(col_index)
                    if is_tuple(col_tup):
                        col_int = map(lambda x:int(x*9+0.49999),col_tup)
                        col_str = string.join(map(lambda x:chr(ord('0')+x),col_int),'')
                        doc_list.append("\\"+col_str+object_name+"\\---")
                    
    key_list = [
        'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10', #,'F11','F12',
        'SHFT-F1','SHFT-F2','SHFT-F3','SHFT-F4','SHFT-F5','SHFT-F6','SHFT-F7',
        'SHFT-F8','SHFT-F9','SHFT-F10']# ,'SHFT-F11','SHFT-F12']
    doc_list.append("")
    doc_list.append("\\888Toggles:")
    zoom_context = 1
                                  
    global labels,ligands,waters,sites,cgos,zooms,dashes
    labels = 1
    ligands = 1
    waters = 1
    sites = 1
    cgos = 0
    zooms = 0
    dashes = 1
    global m4x_sites,m4x_ligands,m4x_waters
    m4x_sites = "m4x_sites"
    m4x_ligands = "m4x_ligands"
    m4x_waters = "m4x_waters"

    cmd.select(m4x_sites,"none")
    cmd.select(m4x_ligands,"none")
    cmd.select(m4x_waters,"none")   
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_zooms)
        doc_list.append(key+": Zoom")
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_sites)
        doc_list.append(key+": Sites")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_waters)
        doc_list.append(key+": Waters")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_dashes)
        doc_list.append(key+": H-Bonds")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_cgos)
        doc_list.append(key+": Fits")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_ligands)
        doc_list.append(key+": Ligands")        
    if len(key_list):
        key = key_list.pop(0)
        cmd.set_key(key,toggle_labels)
        doc_list.append(key+": HB-Dists")        
    
    for a in list:
        include_flag = 0
        water = a+"_water"
        ligand = a+"_ligand"
        site = a+"_site"
        hbond = a+"_hbond"
        if cmd.count_atoms(site):
            if cmd.count_atoms(site+" & m4x_aligned"):
                include_flag = 1
        if cmd.count_atoms(ligand):
            if cmd.count_atoms(ligand+" & m4x_nearby"):
                include_flag = 1
        if include_flag:
            name_list = dict[a]
            if water in name_list:
                cmd.select(m4x_waters,m4x_waters+"|"+water)
            if ligand in name_list:
                cmd.select(m4x_ligands,m4x_ligands+"|"+ligand)
            if site in name_list:
                cmd.select(m4x_sites,m4x_sites+"|"+site+
                "|((byres (neighbor ("+site+" and name c))) and name n+ca)"+
                "|((byres (neighbor ("+site+" and name n))) and name c+ca+o)")
    cmd.wizard("fedora",doc_list)
    toggle_cgos(1)
    toggle_labels(0)
    toggle_dashes(0)
    toggle_ligands(1)
    toggle_sites(0)
    toggle_waters(0)
    toggle_cgos(1)
    cmd.deselect()
#    cmd.feedback("enable","python","output")
    cmd.feedback("enable","objectmolecule","results")
    cmd.set("internal_feedback",1)
    cmd.set("internal_prompt",0)
    cmd.feedback("disable","selector","actions")
    cmd.feedback("disable","scene","actions")
