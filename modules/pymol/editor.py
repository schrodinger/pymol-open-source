import pymol
import cmd

tmp_editor = "_tmp_editor"
tmp_ed_save = "_tmp_ed_save"

# routines to assist in molecular editing

def attach_fragment(selection,fragment,hydrogen,anchor):
   if not selection in cmd.get_names("selections"):
      cmd.fragment(fragment)
   else:
      cmd.fragment(fragment,tmp_editor)
      if cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
         cmd.fuse("(%s and id %d)"%(tmp_editor,hydrogen),"(pk1)",1)
      else:
         cmd.remove("(%s and id %d)"%(tmp_editor,hydrogen))
         cmd.fuse("(%s and id %d)"%(tmp_editor,anchor),"(pk1)",1)
      cmd.delete(tmp_editor)
      
def attach_amino_acid(selection,amino_acid):
   if not selection in cmd.get_names("selections"):
      cmd.fragment(amino_acid)
      if cmd.count_atoms("((%s) and name c)"%amino_acid,quiet=1):
         cmd.edit("((%s) and name c)"%amino_acid)
   else:
      cmd.fragment(amino_acid,tmp_editor)
      if cmd.count_atoms("((%s) and elem n)"%selection,quiet=1):
         cmd.select(tmp_ed_save,"(%s)"%selection)
         cmd.iterate("(%s)"%selection,"stored.resv=resv")
         pymol.stored.resi = str(pymol.stored.resv-1)
         cmd.alter(tmp_editor,"resi=stored.resi")
         cmd.fuse("(%s and name C)"%(tmp_editor),"(pk1)",2)
         sele = ("(name N and (byres neighbor %s) and not (byres %s))"%
                 (tmp_ed_save,tmp_ed_save))
         if cmd.count_atoms(sele,quiet=1):
            cmd.edit(sele)
         cmd.delete(tmp_ed_save)                  
      elif cmd.count_atoms("((%s) and elem c)"%selection,quiet=1):
         cmd.select(tmp_ed_save,"(%s)"%selection)
         cmd.iterate("(%s)"%selection,"stored.resv=resv")
         pymol.stored.resi = str(pymol.stored.resv+1)
         cmd.alter(tmp_editor,"resi=stored.resi")
         cmd.fuse("(%s and name N)"%(tmp_editor),"(pk1)",2)
         sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                 (tmp_ed_save,tmp_ed_save))
         if cmd.count_atoms(sele,quiet=1):
            cmd.edit(sele)
         cmd.delete(tmp_ed_save)         
   cmd.delete(tmp_editor)
      
