
import chempy

class Base:

#------------------------------------------------------------------------------
   def update_index(self):
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
         b = self.bond[i]
         if b.index[0] == index or b.index[1] == index:
            templist.append(i)
      for i in range(len(templist)):
         i = templist[i] - i
         del self.bond[i]

# re-index bond table
      for b in self.bond:
         if b.index[0] > index: 
            b.index[0] = b.index[0] - 1
         if b.index[1] > index: 
            b.index[1] = b.index[1] - 1

#------------------------------------------------------------------------------
   def insert_atom(self,index,atom):

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
   def add_atom(self,atom):
      self.atom.append(atom)
      index = self.nAtom - 1
      if self.index:
         self.index[id(atom)] = index
      return index

#------------------------------------------------------------------------------
   def add_bond(self,bond):
      self.bond.append(bond)      

#------------------------------------------------------------------------------
   def remove_bond(self,index):
      nBond=len(self.Bond)
      del self.bond[index]

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
   def __init__(self,indexed):
      if indexed:
         self.molecule=indexed.molecule
         self.atom = indexed.atom
         self.bond = []
         self.index = None
         for a in self.atom:
            self.bond.append([])
         for b in indexed.bond:
            self.bond[b.index[0]].append(b) # note two refs to same object
            self.bond[b.index[1]].append(b) # note two refs to same object 
         indexed.reset()
      else:
         self.reset()
  
#------------------------------------------------------------------------------
   def reset(self):
      self.index = None
      self.molecule = chempy.Molecule()
      self.atom = []
      self.bond = []
      
#------------------------------------------------------------------------------
   def convert_to_indexed(self):
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
            if b.index[0] == index or b.index[1] == index:
               templist.append(i)
            i = i + 1
         for i in range(len(templist)):
            i = templist[i] - i
            del a[i]

# re-index bond table
      for b in self.bond:
         if b.index[0] > index: 
            b.index[0] = b.index[0] - 1
         if b.index[1] > index: 
            b.index[1] = b.index[1] - 1

#------------------------------------------------------------------------------
   def add_atom(self,atom):
      self.atom.append(atom)
      self.bond.append([])
      index = self.nAtom - 1
      if self.index:
         self.index[id(atom)] = index
      return index

#------------------------------------------------------------------------------
   def sort(self):
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
      



