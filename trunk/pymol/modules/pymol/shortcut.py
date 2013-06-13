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
    try:
        from pymol.checking import is_string, is_list
    except:
        from checking import is_string, is_list
    
    def mkabbr(a, m=1):
        b = a.split('_')
        b[:-1] = [c[0:m] for c in b[:-1]]
        return '_'.join(b)

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
            abbr_dict = self.abbr_dict
            for b in range(1,len(a)):
                sub = a[0:b]
                hash[sub] = 0 if sub in hash else a
            if '_' in a:
                for n in (1, 2):
                    abbr = mkabbr(a, n)
                    if a!=abbr:
                        if abbr in abbr_dict:
                            if a not in abbr_dict[abbr]:
                                abbr_dict[abbr].append(a)
                        else:
                            abbr_dict[abbr]=[a]
                        for b in range(abbr.find('_')+1,len(abbr)):
                            sub = abbr[0:b]
                            hash[sub] = 0 if sub in hash else a

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
            #
            for a in self.keywords:
                for b in range(1,len(a)):
                    sub = a[0:b]
                    hash[sub] = 0 if sub in hash else a
                if '_' in a:
                    for n in (1, 2):
                        abbr = mkabbr(a, n)
                        if a!=abbr:
                            if abbr in abbr_dict:
                                abbr_dict[abbr].append(a)
                            else:
                                abbr_dict[abbr]=[a]
                            for b in range(abbr.find('_')+1,len(abbr)):
                                sub = abbr[0:b]
                                hash[sub] = 0 if sub in hash else a

            self._rebuild_finalize()

        def _rebuild_finalize(self):
            hash = self.shortcut
            for a, adk in self.abbr_dict.iteritems():
                if len(adk)==1:
                    hash[a]=adk[0]
            for a in self.keywords:
                hash[a]=a

        def interpret(self,kee, mode=0):
            '''
            Returns None (no hit), str (one hit) or list (multiple hits)

            kee = str: query string, setting prefix or shortcut
            mode = 0/1: if mode=1, do prefix search even if kee has exact match
            '''
            if not len(kee): # empty string matches everything
                return copy.deepcopy(self.keywords)

            try:
                r = self.shortcut[kee]
            except KeyError:
                return None
            if r and not mode:
                return r

            # prefix search
            lst_set = set(a for a in self.keywords if a.startswith(kee))
            for abbr, a_list in self.abbr_dict.iteritems():
                if abbr.startswith(kee):
                    lst_set.update(a_list)

            # no match
            if not lst_set:
                return None

            # single match: str
            lst = list(lst_set)
            if len(lst) == 1:
                return lst[0]

            # multiple matches: list
            return lst

        def has_key(self,kee):
            return self.shortcut.has_key(kee)

        def __getitem__(self,kee):
            return self.shortcut.get(kee, None)

        def __delitem__(self,kee):
            self.keywords.remove(kee)
            self.rebuild()

        def append(self,kee):
            self.keywords.append(kee)
            self.add_one(kee)
            self._rebuild_finalize()

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
