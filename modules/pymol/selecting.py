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
          _feedback,fb_module,fb_mask, \
          DEFAULT_ERROR, DEFAULT_SUCCESS, _raising, is_ok, is_error
    
    import pymol

    def deselect():
        '''
DESCRIPTION

    "deselect" disables any and all visible selections

USAGE

    deselect

PYMOL API

    cmd.deselect()
        '''
        r = DEFAULT_SUCCESS
        arg = cmd.get_names("selections",enabled_only=1)
        for a in arg:
            cmd.disable(a)
        if _raising(r): raise pymol.CmdException                  
        return r
    

    def select(name,selection="",enable=-1,quiet=1,merge=0):
        '''
DESCRIPTION

    "select" creates a named selection from an atom selection.

USAGE

    select (selection)
    select name, (selection)
    select name = (selection)            # (DEPRECATED)

PYMOL API

    cmd.select(string name, string selection)

EXAMPLES 

    select near , (ll expand 8)
    select near , (ll expand 8)
    select bb, (name ca,n,c,o )

NOTES

    'help selections' for more information about selections.
        '''
        r = DEFAULT_ERROR
        try:
            lock()
            if selection=="":
                selection = name                    
                if _cmd.get("auto_number_selections")!=0.0:
                    sel_cnt = _cmd.get("sel_counter") + 1.0
                    _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
                    name = "sel%02.0f" % sel_cnt
                else:
                    name = "sele"
            else:
                name = name
            # preprocess selection (note: inside TRY)
            selection = selector.process(selection)
            merge = int(merge)
            if merge==1:
                selection = "("+selection+") or ?"+name # merge if exists
            elif merge==2:
                selection = "("+selection+") or ??"+name # merge if exists and active
            #
            r = _cmd.select(str(name),str(selection),int(quiet))
            enable = int(enable)
            if is_ok(r) and enable>0:
                r = _cmd.onoff(str(name),1);
            elif enable == 0:
                r = _cmd.onoff(str(name),0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException                  
        return r


    def pop(name,source,enable=-1,quiet=1):
        r = DEFAULT_ERROR
        try:
            lock()
            r = _cmd.pop(str(name),str(source),int(quiet))
            if is_ok(r):
                enable = int(enable)
                if enable>0:
                    r = _cmd.onoff(str(name),1);
                elif enable == 0:
                    r = _cmd.onoff(str(name),0)
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException                  
        return r

    id_type_dict = {
        'index' : 0,
        'id'    : 1,
        'rank'  : 2,
        }
    
    id_type_sc = Shortcut(id_type_dict.keys())
    
    def select_list(name,object,id_list,state=0,mode='id',quiet=1,):
        '''
DESCRIPTION
    "select_list" is currently in development
    
        '''
        #
        r = DEFAULT_ERROR
        mode = id_type_dict[id_type_sc.auto_err(mode,'identifier type')]
        try:
            lock()
            r = _cmd.select_list(str(name),str(object),list(id_list),int(state)-1,int(mode),int(quiet))
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException
        return r

    def indicate(selection="(all)"):
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
            lock()   
            r = _cmd.select("indicate","("+str(selection)+")",1)
            cmd.enable("indicate")
        finally:
            unlock(r)
        if _raising(r): raise pymol.CmdException                  
        return r







