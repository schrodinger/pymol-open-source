#A* -------------------------------------------------------------------
#B* This file contains source code for the PyMOL computer program
#C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific. 
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

if __name__=='pymol.selecting':
    
    import selector

    import cmd

    from cmd import _cmd,lock,unlock,Shortcut, \
          _feedback,fb_module,fb_mask, is_tuple, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    
    import pymol

    def deselect(_self=cmd):
        '''
DESCRIPTION

    "deselect" disables any and all visible selections

USAGE

    deselect

PYMOL API

    cmd.deselect()
        '''
        r = DEFAULT_SUCCESS
        arg = _self.get_names("selections",enabled_only=1)
        for a in arg:
            _self.disable(a)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r
    

    def select(name, selection="", enable=-1, quiet=1, merge=0, state=0, domain='',_self=cmd): 
        '''
DESCRIPTION

    "select" creates a named atom selection from a
    selection-expression.

USAGE

    select name, selection [, enable [, quiet [, merge [, state [, domain ]]]]]

ARGUMENTS

    name = a unique name for the selection

    selection = a selection-expression

NOTES

    If a selection-expression with explicit surrounding parethenses is
    provided as the first argument, then the default selection name
    is used as the name argument.

EXAMPLES 

    select chA, chain A
    select ( resn his )
    select near142, resi 142 around 5

PYMOL API

    cmd.select(string name, string selection)

SEE ALSO

    delete
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            if selection=="":
                selection = name                    
                if _cmd.get(_self._COb,"auto_number_selections")!=0.0:
                    sel_cnt = _cmd.get(_self._COb,"sel_counter") + 1.0
                    _cmd.legacy_set(_self._COb,"sel_counter","%1.0f" % sel_cnt)
                    name = "sel%02.0f" % sel_cnt
                else:
                    name = "sele"
            if name == None:
                sel_cnt = _cmd.get(_self._COb,"sel_counter") + 1.0
                _cmd.legacy_set(_self._COb,"sel_counter","%1.0f" % sel_cnt)
                name = "sel%02.0f" % sel_cnt
                
            # preprocess selection (note: inside TRY)
            selection = selector.process(selection)
            merge = int(merge)
            if merge==1:
                selection = "("+selection+") or ?"+name # merge if exists
            elif merge==2:
                selection = "("+selection+") or ??"+name # merge if exists and active
            #
            r = _cmd.select(_self._COb,str(name),str(selection),int(quiet),int(state)-1,str(domain))
            enable = int(enable)
            if is_ok(r) and enable>0:
                _cmd.onoff(_self._COb,str(name),1,0);
            elif enable == 0:
                _cmd.onoff(_self._COb,str(name),0,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r


    def pop(name, source, enable=-1, quiet=1, _self=cmd):
        '''
DESCRIPTION

    "pop" provides a mechanism of iterating through an atom selection
    atom by atom, where each atom is sequentially assigned to the
    named selection.
    
USAGE

    pop name, source
    
EXAMPLE

    select src, name ca

    python
    while cmd.pop("tmp","src"):
        cmd.zoom("tmp",2, animate=1)
        for a in range(30):
           cmd.refresh()
           time.sleep(0.05)
    python end
    
PYMOL API

    cmd.deselect()
        '''
        r = DEFAULT_ERROR
        try:
            _self.lock(_self)
            r = _cmd.pop(_self._COb,str(name),str(source),int(quiet))
            if is_ok(r):
                enable = int(enable)
                if enable>0:
                    r = _cmd.onoff(_self._COb,str(name),1,0);
                elif enable == 0:
                    r = _cmd.onoff(_self._COb,str(name),0,0)
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r

    id_type_dict = {
        'index' : 0,
        'id'    : 1,
        'rank'  : 2,
        }
    
    id_type_sc = Shortcut(id_type_dict.keys())
    
    def select_list(name,object,id_list,state=0,mode='id',quiet=1,_self=cmd):
        '''
DESCRIPTION
    "select_list" is currently in development
    
        '''
        #
        r = DEFAULT_ERROR
        mode = id_type_dict[id_type_sc.auto_err(mode,'identifier type')]
        try:
            _self.lock(_self)
            r = _cmd.select_list(_self._COb,str(name),str(object),list(id_list),int(state)-1,int(mode),int(quiet))
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException
        return r

    def indicate(selection="(all)",_self=cmd):
        '''
DESCRIPTION

    "indicate" shows a visual representation of an atom selection.

USAGE

    indicate (selection)

PYMOL API

    cmd.count(string selection)

        '''
        r = DEFAULT_ERROR
        # preprocess selection
        selection = selector.process(selection)
        #      
        try:
            _self.lock(_self)   
            r = _cmd.select(_self._COb,"indicate","("+str(selection)+")",1,-1,'')
            if is_error(r):
                _self.delete("indicate")
            else:
                _self.enable("indicate")
        finally:
            _self.unlock(r,_self)
        if _self._raising(r,_self): raise pymol.CmdException                  
        return r







