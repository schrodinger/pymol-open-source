
import copy
import types
import re
import string

abbr_re = re.compile(r"[^\_]*\_")

def is_string(obj):
   return isinstance(obj,types.StringType)

class Shortcut:

   def __call__(self):
      return self
   
   def __init__(self,list,filter_leading_underscore=1):
      if filter_leading_underscore:
         self.keywords=filter(lambda x:x[:1]!='_',list)
      else:
         self.keywords = copy.deepcopy(list)
      self.shortcut = {}
      self.abbr_dict = {}
      self.rebuild()

   def add_one(self,a):
      # optimize symbols
      hash = self.shortcut
      abbr_dict = self.abbr_dict
      find = string.find
      abbr_re_sub = abbr_re.sub
      for b in range(1,len(a)):
         sub = a[0:b]
         if hash.has_key(sub):
            hash[sub]=0
         else:
            hash[sub]=a
      if find(a,"_")>=0:
         abbr = abbr_re_sub(lambda x:x.group(0)[0]+"_",a)
         if a!=abbr:
            abbr_dict[abbr]=a
            for b in range(string.find(abbr,'_')+1,len(abbr)):
               sub = abbr[0:b]
               if hash.has_key(sub):
                  hash[sub]=0
               else:
                  hash[sub]=a
      
   def rebuild(self):
      # optimize symbols
      hash = self.shortcut
      abbr_dict = self.abbr_dict
      find = string.find
      abbr_re_sub = abbr_re.sub
      #
      for a in self.keywords:
         for b in range(1,len(a)):
            sub = a[0:b]
            if hash.has_key(sub):
               hash[sub]=0
            else:
               hash[sub]=a
         if find(a,"_")>=0:
            abbr = abbr_re_sub(lambda x:x.group(0)[0]+"_",a)
            if a!=abbr:
               abbr_dict[abbr]=a
               for b in range(string.find(abbr,'_')+1,len(abbr)):
                  sub = abbr[0:b]
                  if hash.has_key(sub):
                     hash[sub]=0
                  else:
                     hash[sub]=a
                     
      for a in abbr_dict.keys():
         hash[a]=abbr_dict[a]
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
         for a in self.abbr_dict.keys():
            if a[0:lcm] == kee:
               lst.append(self.abbr_dict[a])
         if(len(lst)==1):
            return lst[0]
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
      self.add_one(kee)
      hash = self.shortcut      
      for a in self.abbr_dict.keys():
         hash[a]=self.abbr_dict[a]
      for a in self.keywords:
         hash[a]=a

   def auto_err(self,kee,descrip=None):
      result = None
      if not self.shortcut.has_key(kee):
         if descrip!=None:
            print "Error: unknown %s: '%s'."%(
               descrip,kee)
            raise parsing.QuietException
      else:
         result = self.shortcut[kee]
         if not is_string(result):
            if descrip!=None:
               print "Error: ambiguous %s:"%descrip
               lst = parsing.list_to_str_list(result)
               for a in lst:
                  print a
               raise QuietException
      return result
              
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

import parsing
