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

   from cmd import _cmd,lock,unlock,Shortcut,QuietException
   from cmd import _feedback,fb_module,fb_mask


   def deselect():
      '''
DESCRIPTION

   "deselect" disables any and all visible selections

USAGE

   deselect

PYMOL API

   cmd.deselect()
      '''
      arg = cmd.get_names("selections")
      for a in arg:
         cmd.disable(a)

   def select(name,selection="",show=-1,quiet=1):
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
      try:
         lock()
         if selection=="":
            sel_cnt = _cmd.get("sel_counter") + 1.0
            _cmd.legacy_set("sel_counter","%1.0f" % sel_cnt)
            selection = name
            name = "sel%02.0f" % sel_cnt
         else:
            name = name
         # preprocess selection (note: inside TRY)
         selection = selector.process(selection)
         #
         r = _cmd.select(str(name),str(selection),int(quiet))
         show = int(show)
         if r and show>0:
            _cmd.onoff(str(name),1);
         elif show == 0:
            _cmd.onoff(str(name),0)
      finally:
         unlock()
      return r

   def pop(name,source,show=-1,quiet=1):
      try:
         lock()
         r = _cmd.pop(str(name),str(source),int(quiet))
         if r<0:
            raise QuietException
         show = int(show)
         if r and show>0:
            _cmd.onoff(str(name),1);
         elif show == 0:
            _cmd.onoff(str(name),0)
      finally:
         unlock()
      return r      

   id_type_dict = {
      'index' : 0,
      'id'    : 1,
      'rank'  : 2,
      }
   
   id_type_sc = Shortcut(id_type_dict.keys())
   
   def select_list(name,object,id_list,id_type='index',show=-1,quiet=1):
      '''
DESCRIPTION
   "select_list" is currently in development
   
      '''
      #
      id_type = id_type_dict[id_type_sc.auto_err(id_type,'identifier type')]
      try:
         lock()
         r = _cmd.select_list(str(name),str(object),list(id_list),int(quiet),int(id_type))
         show = int(show)
         if r and show>0:
            r = _cmd.onoff(str(name),1);
         elif show == 0:
            r = _cmd.onoff(str(name),0)
      finally:
         unlock()   
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
      # preprocess selection
      selection = selector.process(selection)
      #      
      try:
         lock()   
         r = _cmd.select("indicate","("+str(selection)+")",1)
         cmd.enable("indicate")
      finally:
         unlock()
      return r






