
import chempy

class Base:

#------------------------------------------------------------------------------
   def update_index(self):
      if chempy.feedback['verbose']:
         print " "+str(self.__class__)+": updating indexes..."
      c = 0
      self.index = {}
      idx = self.index
      for a in self.atom:
         idx[id(a)] = c
         c = c + 1
         
#------------------------------------------------------------------------------
   def get_residues(self):
      list = []
      if self.nAtom:
         last = self.atom[0]
         c = 0
         start = 0
         for a in self.atom:
            if not a.in_same_residue(last):
               list.append((start,c))
               start = c
               last = a
            c = c + 1
         if (c-start>1):
            list.append((start,c))
      return list

#------------------------------------------------------------------------------

   def get_coord_list(self):
      lst = []
      for a in self.atom:
         lst.append(a.coord)
      return lst

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

class Indexed(Base):

   attr_value = {
      'nAtom' : compile('len(self.atom)','Indexed','eval'),
      'nBond' : compile('len(self.bond)','Indexed','eval'),
      }

   def __getattr__(self,attr):
      if Indexed.attr_value.has_key(attr):
         return eval(Indexed.attr_value[attr])
      else:
         raise AttributeError(attr)
      
#------------------------------------------------------------------------------
   def __init__(self):
      self.reset()
      
#------------------------------------------------------------------------------
   def reset(self):
      self.index = None
      self.molecule = chempy.Molecule()
      self.atom = []
      self.bond = []
   
#------------------------------------------------------------------------------
   def merge(self,other): # steals atom objects from 'other' and resets 'other'
      if chempy.feedback['actions']:
         print " "+str(self.__class__)+": merging models..."
      nAtom=self.nAtom
      self.atom.extend(other.atom)
      for b in other.bond:
         b.index[0]=b.index[0]+nAtom
         b.index[1]=b.index[1]+nAtom
         self.bond.append(b)
      other.reset()
      if self.index:
         self.update_index()

#--------------------------------------------------------------------------------
   def delete_atom(self,index):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+": deleting atom %d." % index
   
      nAtom=self.nAtom

# update index if it exists
      if self.index:
         idx = self.index
         for k in idx.keys():
            if idx[k] > index:
               idx[k] = idx[k] - 1
         del idx[id(self.atom[index])]

# delete atom
      del self.atom[index]

# delete bonds associated with this atom

      nBond = len(self.bond)
      templist = []
      for i in range(nBond):
         if index in self.bond[i].index:
            templist.append(i)
      for i in range(len(templist)):
         j = templist[i] - i
         del self.bond[j]

# re-index bond table
      for b in self.bond:
         if b.index[0] > index: 
            b.index[0] = b.index[0] - 1
         if b.index[1] > index: 
            b.index[1] = b.index[1] - 1

#------------------------------------------------------------------------------
   def insert_atom(self,index,atom):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+': inserting atom "%s" before %d.' % (
            atom.name,index)

      nAtom=self.nAtom
      self.atom.insert(index,atom)
      
# re-index bond table
      for b in self.bond:
         if b.index[0] >= index: 
            b.index[0] = b.index[0] + 1
         if b.index[1] >= index: 
            b.index[1] = b.index[1] + 1

# update index if it exists
      if self.index:
         idx = self.index
         for k in idx.keys():
            if idx[k] >= index:
               idx[k] = idx[k] + 1
         idx[id(atom)] = index
#------------------------------------------------------------------------------
   def index_atom(self,atom):
      c = 0
      id_at = id(atom)
      for a in self.atom:
         if id(a)==id_at:
            return c
         c = c + 1
      return -1
      
#------------------------------------------------------------------------------
   def add_atom(self,atom):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+': adding atom "%s".' % atom.name
      self.atom.append(atom)
      index = self.nAtom - 1
      if self.index:
         self.index[id(atom)] = index
      return index

#------------------------------------------------------------------------------
   def add_bond(self,bond):
      if chempy.feedback['bonds']:
         print " "+str(self.__class__)+": adding bond (%d,%d)." % \
               (bond.index[0],bond.index[1])
      self.bond.append(bond)      

#------------------------------------------------------------------------------
   def remove_bond(self,index):
      if chempy.feedback['bonds']:
         print " "+str(self.__class__)+": removing bond %d." % index
      nBond=len(self.Bond)
      del self.bond[index]
#------------------------------------------------------------------------------
   def convert_to_connected(self):
      if chempy.feedback['verbose']:
         print " "+str(self.__class__)+": converting to connected model..."
      model = Connected()
      model.molecule = self.molecule
      model.atom = self.atom
      model.bond = []
      model.index = None
      for a in model.atom:
         model.bond.append([])
      for b in self.bond:
         model.bond[b.index[0]].append(b) # note two refs to same object
         model.bond[b.index[1]].append(b) # note two refs to same object 
      self.reset()
      return model
#------------------------------------------------------------------------------
   def from_molobj(self,molobj): 
      self.reset()
      mol = self.molecule
      if len(molobj.title):
         mol.title = molobj.title
      if len(molobj.comments):
         mol.comments = molobj.comments
      mol.chiral = molobj.chiral
      mol.dim_code = molobj.dimcode
      for a in molobj.atom:
         at = chempy.Atom()
         at.symbol = a.symbol
         at.name = a.name
         if a.resn != chempy.Atom.defaults['resn']:
            at.resn = a.resn
         if a.resn_code != chempy.Atom.defaults['resn_code']:
            at.resn_code = a.resn_code
         at.resi = a.resi
         at.resi_number = a.resi_number
         at.b = a.b
         at.q = a.q
         at.alt = a.alt
         at.hetatm = a.hetatm
         if a.segi != chempy.Atom.defaults['segi']:
            at.segi = a.segi
         if a.chain != chempy.Atom.defaults['chain']:
            at.chain = a.chain
         at.color_code = a.color_code
         at.coord = a.coord
         at.formal_charge = a.formal_charge
         at.partial_charge = a.partial_charge
         if a.numeric_type != -99:
            at.numeric_type = a.numeric_type
         if a.text_type != 'UNKNOWN':
            at.text_type = a.text_type
         at.stereo = a.stereo
         if hasattr(a,'flags'):
            at.flags = a.flags
         if hasattr(a,'vdw'):
            at.vdw = a.vdw
         self.atom.append(at)
      for b in molobj.bond:
         bnd = chempy.Bond()
         bnd.index = [b.atom[0],b.atom[1]]
         bnd.order = b.order
         bnd.stereo = b.stereo
         self.bond.append(bnd)
#------------------------------------------------------------------------------
   def sort(self):
      if chempy.feedback['verbose']:
         print " "+__name__+": sorting..."
      if not self.index:
         self.update_index()
      old_index = self.index
      self.atom.sort()      
      self.update_index()
      xref = {}
      new_index = self.index
      for a in new_index.keys():
         xref[old_index[a]] = new_index[a]
      for b in self.bond:
         b.index[0] = xref[b.index[0]]
         b.index[1] = xref[b.index[1]]
      del old_index
      del xref
         
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
class Connected(Base):

   attr_value = {
      'nAtom' : compile('len(self.atom)','Connected','eval'),
      }

   def __getattr__(self,attr):
      if Connected.attr_value.has_key(attr):
         return eval(Connected.attr_value[attr])
      else:
         raise AttributeError(attr)
      
#------------------------------------------------------------------------------
   def __init__(self):
      self.reset()
  
#------------------------------------------------------------------------------
   def reset(self):
      self.index = None
      self.molecule = chempy.Molecule()
      self.atom = []
      self.bond = []
      
#------------------------------------------------------------------------------
   def convert_to_indexed(self):
      if chempy.feedback['verbose']:
         print " "+str(self.__class__)+": converting to indexed model..."
      indexed = Indexed()
      indexed.atom = self.atom
      indexed.molecule = self.molecule
      c = 0 
      for a in self.bond:
         for b in a:
            if b.index[0] == c:
               indexed.bond.append(b)
         c = c + 1
      self.reset()
      return indexed

#------------------------------------------------------------------------------
   def insert_atom(self,index,atom):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+': inserting atom "%s" before %d.' % (
            atom.name,index)

      nAtom=self.nAtom
      self.atom.insert(index,atom)
      
# re-index bond table
      for a in self.bonds:
         for b in a:
            if b.index[0] >= index:
               b.index[0] = b.index[0] + 1
            if b.index[1] >= index:
               b.index[1] = b.index[1] + 1

# update index if it exists
      if self.index:
         idx = self.index
         for k in idx.keys():
            if idx[k] >= index:
               idx[k] = idx[k] + 1
         idx[id(atom)] = index

#------------------------------------------------------------------------------
   def delete_atom(self,index):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+": deleting atom %d." % index

      nAtom=self.nAtom

# update index if it exists
      if self.index:
         idx = self.index
         for k in idx.keys():
            if idx[k] > index:
               idx[k] = idx[k] - 1
         del idx[id(self.atom[index])]

# delete atom
      del self.atom[index]

# delete bonds associated with this atom

      nBond = len(self.bond)

      for a in self.bond:
         i = 0
         templist = []
         for b in a:
            if index in b.index:
               templist.append(i)
            i = i + 1
         for i in range(len(templist)):
            j = templist[i] - i
            del a[j]

# re-index bond table
      for b in self.bond:
         if b.index[0] > index: 
            b.index[0] = b.index[0] - 1
         if b.index[1] > index: 
            b.index[1] = b.index[1] - 1

#------------------------------------------------------------------------------
   def add_atom(self,atom):
      if chempy.feedback['atoms']:
         print " "+str(self.__class__)+': adding atom "%s".' % atom.name
      self.atom.append(atom)
      self.bond.append([])
      index = self.nAtom - 1
      if self.index:
         self.index[id(atom)] = index
      return index

#------------------------------------------------------------------------------
   def sort(self):
      if chempy.feedback['verbose']:
         print " "+__name__+": sorting..."
      if not self.index:
         self.update_index()
      old_index = self.index
      self.atom.sort()      
      self.update_index()
      xref = {}
      new_index = self.index
      for a in new_index.keys():
         xref[old_index[a]] = new_index[a]
      new_bond = [None] * self.nAtom
      c = 0
      tmp_list = []
      for a in self.bond:
         for b in a:
            if c==b.index[0]:
               tmp_list.append(b)
         new_bond[xref[c]] = a
         c = c + 1
      for b in tmp_list:
         b.index[0] = xref[b.index[0]]
         b.index[1] = xref[b.index[1]]
      del self.bond
      self.bond = new_bond
      del old_index
      del xref
      



