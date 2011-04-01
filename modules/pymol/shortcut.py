#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* Copyright (c) Schrodinger, LLC. 
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information. 
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-* 
#-* 
#-*
#Z* -------------------------------------------------------------------

if __name__=='pymol.shortcut' or __name__=='shortcut':
    
    import copy
    import types
    import re
    import string
    
    abbr_re = re.compile(r"[^\_]*\_")

    def is_string(obj):
        return isinstance(obj,types.StringType)

    def is_list(obj):
        return isinstance(obj,types.ListType)

    class Shortcut:

        def __call__(self):
            return self

        def __init__(self,list=[],filter_leading_underscore=1):
            self.filter_leading_underscore = filter_leading_underscore
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
            abbr_dict = {}
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
                    if(abbr_dict.has_key(abbr)):
                        abbr_dict[abbr].append(a)
                    else:
                        abbr_dict[abbr]=[a]
                    for b in range(string.find(abbr,'_')+1,len(abbr)):
                        sub = abbr[0:b]
                        if hash.has_key(sub):
                            hash[sub]=0
                        else:
                            hash[sub]=a

                abbr = abbr_re_sub(lambda x:x.group(0)[0:2]+"_",a)
                if a!=abbr:
                    if(abbr_dict.has_key(abbr)):
                        abbr_dict[abbr].append(a)
                    else:
                        abbr_dict[abbr]=[a]
                    for b in range(string.find(abbr,'_')+1,len(abbr)):
                        sub = abbr[0:b]
                        if hash.has_key(sub):
                            hash[sub]=0
                        else:
                            hash[sub]=a

        def rebuild(self,list=None):
            if list!=None:
                if self.filter_leading_underscore:
                    self.keywords=filter(lambda x:x[:1]!='_',list)
                else:
                    self.keywords = copy.deepcopy(list)
            # optimize symbols
            self.shortcut = {}
            hash = self.shortcut
            self.abbr_dict = {}
            abbr_dict = self.abbr_dict
            self.abbr_list = []
            abbr_list = self.abbr_list
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
                        if(abbr_dict.has_key(abbr)):
                            abbr_dict[abbr].append(a)
                        else:
                            abbr_dict[abbr]=[a]
                        for b in range(string.find(abbr,'_')+1,len(abbr)):
                            sub = abbr[0:b]
                            if hash.has_key(sub):
                                hash[sub]=0
                            else:
                                hash[sub]=a
                    abbr = abbr_re_sub(lambda x:x.group(0)[0:2]+"_",a)
                    if a!=abbr:
                        if(abbr_dict.has_key(abbr)):
                            abbr_dict[abbr].append(a)
                        else:
                            abbr_dict[abbr]=[a]
                        for b in range(string.find(abbr,'_')+1,len(abbr)):
                            sub = abbr[0:b]
                            if hash.has_key(sub):
                                hash[sub]=0
                            else:
                                hash[sub]=a

            for a in abbr_dict.keys():
                adk = abbr_dict[a]
                if len(adk)==1:
                    hash[a]=adk[0]
            for a in self.keywords:
                hash[a]=a

        def interpret(self,kee):
            if not len(kee): # empty string matches everything
                return copy.deepcopy(self.keywords)
            elif not self.shortcut.has_key(kee):
                return None # unrecognized, returns None
            elif self.shortcut[kee]==0:
                lst_dict = {}
                lcm = len(kee)
                for a in self.keywords:
                    if a[0:lcm] == kee:
                        lst_dict[a]=1 # ambiguous returns list
                for a in self.abbr_dict.keys():
                    if a[0:lcm] == kee:
                        for b in self.abbr_dict[a]:
                            lst_dict[b]=1
                lst = lst_dict.keys()
                if(len(lst)==1):
                    return lst[0]
                elif not len(lst):
                    return None
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
                adk = self.abbr_dict[a]
                if len(adk)==1:
                    hash[a]=adk[0]
            for a in self.keywords:
                hash[a]=a

        def auto_err(self,kee,descrip=None):
            result = None
            if not self.shortcut.has_key(kee):
                if descrip!=None:
                    print "Error: unknown %s: '%s'."%(
                        descrip,kee),
                    lst = self.interpret('')
                    if is_list(lst):
                        if len(lst)<100:
                            lst.sort()
                            print "Choices:"
                            lst = parsing.list_to_str_list(lst)
                            for a in lst: print a
                        else:
                            print
                    else:
                        print
                    raise parsing.QuietException

            else:
                result = self.interpret(kee)
                if not is_string(result):
                    if descrip!=None:
                        print "Error: ambiguous %s:"%descrip
                        lst = parsing.list_to_str_list(result)
                        for a in lst:
                            print a
                        raise parsing.QuietException
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
