
import chempy
import copy
from chempy.models import Indexed

from Numeric import *
from Precision import *

class FastModel:

#------------------------------------------------------------------------------
    def __init__(self):
        self.reset()

#------------------------------------------------------------------------------
    def reset(self):
        self.nAtom = 0
        self.molecule = chempy.Molecule()
        self.txta = None
        self.inta = None
        self.flta = None
        self.crda = None
        self.bnda = None

#------------------------------------------------------------------------------
    def from_indexed(self,model):
        self.reset()
        self.nAtom = model.nAtom
        self.nBond = model.nBond
        self.molecule = copy.deepcopy(model.molecule)
        self.txta = resize(array(' ','c'),(self.nAtom,as_width))
        self.inta = zeros((self.nAtom,ai_width),Int32)
        self.flta = zeros((self.nAtom,af_width),Float)
        self.bnda = zeros((self.nBond,bi_width),Int32)
        self.crda = zeros((self.nAtom,3),Float)

        c = 0
        for a in model.atom:
            txt = "%-4s%-2s%-1s%-4s%-1s%-4s%-1s%-4s%-20s" % \
                    (a.name[0:4],a.symbol[0:2],a.alt[0:1],a.resn[0:4],
                     a.resn_code[0:1],a.resi[0:4],a.chain[0:1],
                     a.segi[0:4],a.text_type[0:20])
            self.txta[c] = txt
            self.inta[c] = [ a.resi_number, a.hetatm, a.formal_charge,
                                  a.flags, a.color_code, a.stereo, a.numeric_type ]
            self.flta[c] = [ a.b, a.q, a.partial_charge, a.vdw ]
            self.crda[c] = [ a.coord[0], a.coord[1], a.coord[2] ]
            c = c + 1

        c = 0
        for b in model.bond:
            self.bnda[c] = [ b.index[0], b.index[1], b.order, b.stereo ]
            c = c + 1

#------------------------------------------------------------------------------
    def convert_to_indexed(self):
        
        model = Indexed()
        model.molecule = copy.deepcopy(self.molecule)
        
        for c in xrange(self.nAtom):
            at = chempy.Atom()
            txta = self.txta[c]

            for attrib in ( 'name', 'symbol', 'resn', 'resn_code', 'resi',
                                 'alt', 'chain', 'segi', 'text_type' ):
                ll = as[attrib]
                setattr(at,attrib,string.strip(string.join(txta[ll[0]:ll[1]],'')))

            inta = self.inta[c]
            for attrib in ( 'resi_number', 'hetatm', 'formal_charge','flags',
                                 'color_code', 'stereo', 'numeric_type' ):
                setattr(at,attrib,inta[ai[attrib]])

            flta = self.flta[c]
            for attrib in ( 'b', 'q', 'partial_charge' ) :
                setattr(at,attrib,flta[af[attrib]])

            crda = self.crda[c]
            at.coord = [crda[0],crda[1],crda[2]]

            # probably need to add some checking here to eliminate values
            # which come back as defaults
            
            model.atom.append(at)
            
        for c in xrange(self.nBond):
            bnd = chempy.Bond()
            bnda = self.bnda[c]
            bnd.index = [bnda[bi_index0],bnda[bi_index1]]
            bnd.order = bnda[bi_order]
            bnd.stereo = bnda[bi_stereo]
            model.bond.append(bnd)

        return model
    
#------------------------------------------------------------------------------

# text properties

as = {}
as_width = 0

as['name'] = [ as_width ]
as_width = as_width + 4
as['name'].append(as_width)

as['symbol'] = [ as_width ]
as_width = as_width + 2
as['symbol'].append(as_width)

as['alt'] = [ as_width ]
as_width = as_width + 1
as['alt'].append(as_width)

as['resn'] = [ as_width ]
as_width = as_width + 4
as['resn'].append(as_width)

as['resn_code'] = [ as_width ]
as_width = as_width + 1
as['resn_code'].append(as_width)

as['resi'] = [ as_width ]
as_width = as_width + 4
as['resi'].append(as_width)

as['chain'] = [ as_width ]
as_width = as_width + 1
as['chain'].append(as_width)

as['segi' ] = [ as_width ]
as_width = as_width + 4
as['segi'].append(as_width)

as['text_type'] = [ as_width ]
as_width = as_width + 20
as['text_type'].append(as_width)

# integer properties

ai = {}
ai_width = 0
ai['resi_number']    = ai_width
ai_width = ai_width + 1
ai['hetatm']         = ai_width
ai_width = ai_width + 1
ai['formal_charge']  = ai_width
ai_width = ai_width + 1
ai['flags']          = ai_width
ai_width = ai_width + 1
ai['color_code']     = ai_width
ai_width = ai_width + 1
ai['stereo']         = ai_width
ai_width = ai_width + 1
ai['numeric_type']   = ai_width
ai_width = ai_width + 1

# float properties

af = {}
af_width = 0
af['b']              = af_width
af_width = af_width + 1
af['q']              = af_width
af_width = af_width + 1
af['partial_charge'] = af_width
af_width = af_width + 1
af['vdw'] = af_width
af_width = af_width + 1

# bond information

bi_index0         = 0
bi_index1         = 1
bi_order          = 2
bi_stereo         = 3
bi_width          = 4


