try:
   import _champ
except ImportError:
   print " Error: unable to import architecture-dependent _champ module."

# okay, this module is going to take some planning, since Champ
# has the potential to evolve into a full-blown
# chemical informatics system 

# what do we need to be prepared to handle?
#
# - single structures
# - single patterns
# - lists of patterns
# - lists of compounds
# - dictionaries of compounds (key/value)
# - molecules in other class formats (chempy models, etc.)
# - subsets of protein structures
# - scripts which apply rules according to patterns
# - synchronized chempy models

# what does the architecture-dependent modules need to support?

# - Smiles/SMARTS/etc. conversion (bidirectional)
# - ChemPy model conversion (bidirectional)
# - modification recipe for extant chempy models?...

class Champ:

   def __init__(self):
      '''
      creates a new champ object, which contains
      a collect of comparable, searchable patterns 
      '''
      self._champ = _champ.new()

   def insert_pattern_string(self,pattern):
      '''
      adds a new pattern, and
      returns the index to this pattern
      '''
      (e,r) = _champ.insert_pattern_string(self._champ,str(pattern))
      if e: raise RuntimeError
      return r

   def pattern_free(self,index):
      '''
      frees a pattern
      '''
      (e,r) = _champ.pattern_free(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_clear_tags(self,index):
      '''
      resets a pattern\'s tags
      '''
      (e,r) = _champ.pattern_clear_tags(self._champ,int(index))
      if e: raise RuntimeError
      return r
      
   def get_pattern_string(self,index):
      '''
      retrieves the smiles string for a given pattern index
      '''
      (e,r) = _champ.get_pattern_string(self._champ,int(index))
      if e: raise RuntimeError
      return r
   
   def insert_model(self,model):
      '''
      inserts the pattern from a ChemPy model pattern and returns the index
      '''
      (e,r) = _champ.insert_model(self._champ,model)
      if e: raise RuntimeError
      return r

   def list_prepend_pattern_strings(self,handle,pattern):
      '''
      adds pattern string to the top of a list
      and returns the new list identifier
      '''
      (e,r) = _champ.list_prepend_pattern_strings(self._champ,int(handle),pattern)
      if e: raise RuntimeError
      return r

   def list_new(self):
      '''
      returns a new list handle
      '''
      (e,r) = _champ.list_new(self._champ)
      if e: raise RuntimeError
      return r

   def list_free(self,handle,purge=1):
      '''
      destroys a list and associated patterns (if purge = 1)
      '''
      (e,r) = _champ.list_free(self._champ,int(handle),int(purge))
      if e: raise RuntimeError
      return r

   def list_get_pattern_indices(self,handle):
      '''
      returns pattern indices in a list as a list
      '''
      (e,r) = _champ.list_get_pattern_list(self._champ,int(handle))
      if e: raise RuntimeError
      return r

   def list_get_pattern_strings(self,handle):
      '''
      returns list of smiles strings contained in a list
      '''
      (e,r) = _champ.list_get_pattern_strings(self._champ,int(handle))
      if e: raise RuntimeError
      return r

   def match_1v1_b(self,pattern,target): # boolean
      '''
      returns whether or not pattern can be found in target
      as a boolean result
      '''
      (e,r) = _champ.match_1v1_b(self._champ,
                           int(pattern),int(target))
      if e: raise RuntimeError
      return r

   def match_1v1_n(self,pattern,target,limit,tag=0):
      '''
      returns tag list for single pattern on target
      '''
      (e,r) = _champ.match_1v1_n(self._champ,
                           int(pattern),int(target),
                                 int(limit),int(tag))
      if e: raise RuntimeError
      return r

   def match_1v1_map(self,pattern,target,limit): # boolean
      '''
      returns mappings (if any) between two patterns
      '''
      (e,r) = _champ.match_1v1_map(self._champ,
                           int(pattern),int(target),int(limit))
      if e: raise RuntimeError
      return r

   def match_Nv1_n(self,pattern,target,limit,tag):
      '''
      returns tag list for list of patterns onto a target
      '''
      (e,r) = _champ.match_Nv1_n(self._champ,
                           int(pattern),int(target),
                                 int(limit),int(tag))
      if e: raise RuntimeError
      return r

      
   def match_1vN_n(self,pattern,handle):
      '''
      returns count of how many times pattern occurs in list
      '''
      (e,r) = _champ.match_1vN_n(self._champ,
                           int(pattern),int(handle))
      if e: raise RuntimeError
      return r

   def memory_dump(self):
      '''
      dump bulk memory information
      '''
      (e,r) = _champ.memory_dump(self._champ)
      if e: raise RuntimeError
      return r

   def pattern_get_cycle(self,index):
      '''
      debugging routine
      '''
      (e,r) = _champ.pattern_get_cycle(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_get_class(self,index):
      '''
      debugging routine
      '''
      (e,r) = _champ.pattern_get_class(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_get_codes(self,index):
      '''
      debugging routine
      '''
      (e,r) = _champ.pattern_get_codes(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_get_tags(self,index):
      '''
      get atomic tags
      '''
      (e,r) = _champ.pattern_get_tags(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_get_atom_symbols(self,index):
      '''
      debugging routine
      '''
      (e,r) = _champ.pattern_get_atom_symbols(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_dump(self,index):
      '''
      debugging routine
      '''
      (e,r) = _champ.pattern_dump(self._champ,int(index))
      if e: raise RuntimeError
      return r

   def pattern_orient_bonds(self,index):
      '''
      experimental
      '''
      (e,r) = _champ.pattern_orient_bonds(self._champ,int(index))
      if e: raise RuntimeError
      return r
      
