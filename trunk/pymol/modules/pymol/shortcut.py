
import copy

class Shortcut:
   
   def __init__(self,list):
      self.keywords = copy.deepcopy(list)
      self.shortcut = {}
      self.rebuild()

   def rebuild(self):
      hash = self.shortcut
      for a in self.keywords:
         for b in range(1,len(a)):
            sub = a[0:b]
            if hash.has_key(sub):
               hash[sub]=0
            else:
               hash[sub]=a
      for a in self.keywords:
         hash[a]=a

   def interpret(self,kee):
      if not len(kee): # empty string matches everything
         return copy.deepcopy(self.keywords)
      elif not self.shortcut.has_key(kee):
         return None # unrecognized, returns None
      elif self.shortcut[kee]==0:
         lst = []
         lcm = len(kee)
         for a in self.keywords:
            if a[0:lcm] == kee:
               lst.append(a) # ambiguous returns list
         return lst
      else:
         return self.shortcut[kee] # otherwise return string
      
   def has_key(self,kee):
      return self.shortcut.has_key(kee)

   def __getitem__(self,kee):
      if self.shortcut.has_key(kee):
         return self.shortcut[kee]
      else:
         return None

   def __delitem__(self,kee):
      self.keywords.remove(kee)
      self.rebuild()
      
   def append(self,kee):
      self.keywords.append(kee)
      hash = self.shortcut
      for b in range(1,len(kee)+1):
         sub = a[0:b]
         if hash.has_key(sub):
            hash[sub]=0
         else:
            hash[sub]=kee
      for a in self.keywords:
         hash[a]=a

if __name__=='__main__':
   list = [ 'warren','wasteland','electric','well' ]
   sc = Shortcut(list)
   tv = sc.has_key('a')
   print tv==0,tv
   tv = sc.has_key('w')
   print tv==1,tv
   tv = sc.has_key('war')
   print tv==1,tv
   tv = sc.ambiguous('w')
   print tv==['warren', 'wasteland', 'well'],tv   
   tv = sc.ambiguous('e')
   print tv==None,tv   
