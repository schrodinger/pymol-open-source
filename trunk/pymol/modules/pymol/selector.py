import re
import string
import traceback

pat_re = re.compile(r'[^ \(\)\!\&\|]*\/[^ \(\)\!\&\|]*')
num_re = re.compile(r'1|2|3|4|5|6|7|8|9|0')

def process(sele): # expand slash notation into a standard atom selection

   sele = str(sele)
   if string.find(sele,'/')<0:
      return sele
   while 1:
      mo = pat_re.search(sele)
      if mo == None:
         break
      bef = sele[:mo.start()]
      aft = sele[mo.end():]
      mat = mo.group()
      arg = string.split(mat,"/")
      la = len(arg)
      # wildcard -> blank
      arg = map(string.strip,arg)
      for c in range(len(arg)):
         if arg[c]=='*': arg[c]=''
      # defaults
      model = ''
      segment = ''
      chain = ''
      residue = ''
      name = ''
      # interpret
      if mat[0]=="/": # preceeding slash, interpret left to right
         if la>1: model = arg[1]
         if la>2: segment = arg[2]
         if la>3: chain = arg[3]
         if la>4: residue = arg[4]
         if la>5: name = arg[5]
      else: # no preceeding slash, interpret right to left
         arg.reverse()
         if la>0: name = arg[0]
         if la>1: residue = arg[1]
         if la>2: chain = arg[2]
         if la>3: segment = arg[3]
         if la>4: model = arg[4]
      # lst
      lst = []
      if model!='': lst.append("("+string.replace(model,"+","|")+")")
      if segment!='': lst.append("s;"+string.replace(segment,'+',','))
      if chain!='': lst.append("c;"+string.replace(chain,'+',','))
      if residue!='':
         if num_re.search(residue)==None:
            lst.append("n;"+residue)            
         else:
            residue = string.replace(residue,'+',',')
            residue = string.replace(residue,'-',':')
            if ((string.find(residue,',')>=0) and # compound residue specification
                (string.find(residue,':')>=0)):
               new_list = []
               for a in string.split(residue,','): # spread it out...
                  new_list.append("i;"+a)
               residue = "("+string.join(new_list,'|')+")"
               lst.append(residue)
            else:
               lst.append("i;"+residue)                        
      if name!='': lst.append("n;"+string.replace(name,'+',','))
      if not len(lst):
         st = "*"
      else:
         st = reduce(lambda x,y:x+"&"+y,lst)
      sele = bef+"("+st+")"+aft
   return sele

