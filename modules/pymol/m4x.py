# 
# Author: Scott Dixon, Metaphorics, LLC
# This source code is contributed to the public domain and may be freely
# copied and distributed for research, profit, fun or any other reason,
# with these restrictions: (1) unmodified or functionally equivalent code
# derived from this code must contain this notice, (2) all derived code
# must acknowledge the author and institution, and (3) no liability is
# assumed by the author(s) for any use or misuse of this software.

import cmd
from cmd import _cmd

from chempy import Storage,Atom,Bond
from chempy.models import Indexed

import string
from chempy import cex as CEX

class CEXpyParser(CEX.CEXsmilesParser):
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
def readcex(file,*args):  
    import os.path

      
    standardResidues = ("ALA", "ARG", "ASP", "ASN", "ASX", "CYS", "GLY", "GLU",
                        "GLN", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
                        "SER", "THR", "TYR", "TRP", "VAL", "HID", "HIE")
    def childprop(tag, tree):
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
    _cmd.finish_object(modelname)
    f.close()

def colorbyB(selection="spheres",first=7,last=3):
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
