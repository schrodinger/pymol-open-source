import re
import string
import traceback
import types
        
sele_re = re.compile(r'^[^\(].*[ \(\)\!\&\|]') # unwrapped selection?
pat_re = re.compile(r'[^ \(\)\!\&\|]*\/[^ \(\)\!\&\|]*')
num_re = re.compile(r'1|2|3|4|5|6|7|8|9|0')
obj_name_re = re.compile(r'([A-Za-z0-9_\+\-]+)')

def is_tuple(obj):
    return isinstance(obj,types.TupleType)

def work_around_model(model): # ugly workaround to support use of 'model' as object name
    list = obj_name_re.split(model)
    new_list = []
    for entry in list:
        if entry == 'model':
            new_list.append("(object model)")
        else:
            new_list.append(entry)
    return string.join(new_list,'')
    
def process(sele): # expand slash notation into a standard atom selection
    # convert object/index tuples into selection strings
    if is_tuple(sele):
        sele="%s`%d"%sele
    # convert unicode hyphens to dashes
    sele = string.replace(str(sele),u'\u2212','-')
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
        alt = ''
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
        if model!='':
            model = (model.replace("+","|")).replace("\\|","+")
            model = work_around_model(model)
            lst.append("("+model+")")
        if segment!='': lst.append("s;"+string.replace(segment,'+',','))
        if chain!='': lst.append("c;"+string.replace(chain,'+',','))
        if residue!='':
            res_name = ''
            res_id = ''
            res_split=string.find(residue,'`')
            if res_split>=0: # are resn or resi explicitly indicated using a forward apostrophe?
                res_name = residue[0:res_split]
                res_id = residue[res_split+1:]
            elif num_re.search(residue)==None: # if residue has no numeric character, then treat as residue name...
                res_name = residue
            else: # otherwise treat as a residue identifier...
                res_id = residue
            if len(res_name):
                lst.append("r;"+string.replace(res_name,'+',','))            
            if len(res_id):
                lst.append("i;"+res_id)
#               res_id = string.replace(res_id,'+',',')
#               res_id = string.replace(res_id,'-',':')
#               if ((string.find(res_id,',')>=0) and # compound residue specification
#                   (string.find(res_id,':')>=0)):
#                  new_list = []
#                  for a in string.split(res_id,','): # spread it out...
#                     new_list.append("i;"+a)
#                  res_id = "("+string.join(new_list,'|')+")"
#                  lst.append(res_id)
#               else:
#                  lst.append("i;"+res_id)                        
        if name!='':
            if(string.find(name,'`')>=0): # alternate conformations present
                (name,alt) = string.split(name,'`')
                if not len(alt):
                    alt="''"
                if len(name):
                    lst.append("n;"+name)
            else:
                lst.append("n;"+name)
        if alt!='': lst.append("alt "+string.replace(alt,'+',','))      
        if not len(lst):
            st = "*"
        else:
            st = reduce(lambda x,y:x+"&"+y,lst)
        sele = bef+"("+st+")"+aft
    return sele

