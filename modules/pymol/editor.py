import setting
import pymol
import cmd
import parsing

QuietException = parsing.QuietException

tmp_editor = "_tmp_editor"
tmp_ed_save = "_tmp_ed_save"
tpk1 = "_tmp_tpk1"
tpk2 = "_tmp_tpk2"

# routines to assist in molecular editing

def attach_fragment(selection,fragment,hydrogen,anchor):
   if not selection in cmd.get_names("selections"):
      if fragment in cmd.get_names("objects"):
         print " Error: an object with than name already exists"
         raise QuietException
      else:
         cmd.fragment(fragment)
         if cmd.get_setting_legacy("auto_remove_hydrogens"):
            cmd.remove("(hydro and %s)"%fragment)
   else:
      cmd.fragment(fragment,tmp_editor)
      if cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
         cmd.fuse("(%s and id %d)"%(tmp_editor,hydrogen),"(pk1)",1)
         if cmd.get_setting("auto_remove_hydrogens"):
            cmd.remove("(hydro and pkchain)")            
      else:
         cmd.remove("(%s and id %d)"%(tmp_editor,hydrogen))
         cmd.fuse("(%s and id %d)"%(tmp_editor,anchor),"(pk1)",1)
         if cmd.get_setting("auto_remove_hydrogens"):
            cmd.remove("(hydro and pkchain)")            
      cmd.delete(tmp_editor)
      
def attach_amino_acid(selection,amino_acid):
   if not selection in cmd.get_names("selections"):
      if amino_acid in cmd.get_names("objects"):
         print " Error: an object with than name already exists"
         raise QuietException
      cmd.fragment(amino_acid)
      if cmd.get_setting_legacy("auto_remove_hydrogens"):
         cmd.remove("(hydro and %s)"%amino_acid)
      if cmd.count_atoms("((%s) and name c)"%amino_acid,quiet=1):
         cmd.edit("((%s) and name c)"%amino_acid)
   else:
      ss = int(cmd.get_setting("secondary_structure"))
      if ss:
         if ss==1: # helix
            phi=-57.0
            psi=-47.0
         elif ss==2: # antipara-beta
            phi=-139.0
            psi=135.0
         elif ss==3: # para-beta
            phi=-119.0
            psi=113.0
         else:
            phi=180.0
            psi=180.0

      cmd.fragment(amino_acid,tmp_editor)
      if cmd.count_atoms("((%s) and elem n)"%selection,quiet=1):
         cmd.select(tmp_ed_save,"(%s)"%selection)
         cmd.iterate("(%s)"%selection,"stored.resv=resv")
         pymol.stored.resi = str(pymol.stored.resv-1)
         cmd.alter(tmp_editor,"resi=stored.resi")
         cmd.fuse("(%s and name C)"%(tmp_editor),"(pk1)",2)
         if cmd.get_setting("auto_remove_hydrogens"):
            cmd.remove("(pkchain and hydro)")
         cmd.set_dihedral("(name ca and neighbor pk2)","(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
         cmd.set_geometry("pk2",3,3) # make nitrogen planer
         if ss:
            cmd.select(tpk1,"pk2")
            cmd.select(tpk2,"pk1")
            if amino_acid[0:3]!='pro':
               cmd.set_dihedral( # PHI
                  "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
                  "(name ca and neighbor "+tpk1+")", # CA 
                  tpk1, # N
                  tpk2, # C
                  phi)
            cmd.set_dihedral( # PSI (n-1)
               tpk1, # N
               tpk2, # C
               "(name ca and neighbor "+tpk2+")", # CA
               "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
               psi)
            cmd.delete(tpk1)
            cmd.delete(tpk2)               
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
         if cmd.get_setting("auto_remove_hydrogens"):
            cmd.remove("(pkchain and hydro)")
         cmd.set_dihedral("(name ca and neighbor pk2)","(pk2)","(pk1)","(name ca,ch3 and neighbor pk1)",180.0)
         cmd.set_geometry("pk1",3,3) # make nitrogen planer
         if ss:
            cmd.select(tpk1,"pk1")
            cmd.select(tpk2,"pk2")
            if amino_acid[0:3]!='pro':
               cmd.set_dihedral( # PHI
                  tpk2, # C
                  tpk1, # N
                  "(name ca and neighbor "+tpk1+")", # CA 
                  "(name c and neighbor (name ca and neighbor "+tpk1+"))", # C
                  phi)
            cmd.set_dihedral( # PSI (n-1)
               "(name n and neighbor (name ca and neighbor "+tpk2+"))", # C
               "(name ca and neighbor "+tpk2+")", # CA
               tpk2, # C
               tpk1, # N
               psi)
            cmd.delete(tpk1)
            cmd.delete(tpk2)               
         sele = ("(name C and (byres neighbor %s) and not (byres %s))"%
                 (tmp_ed_save,tmp_ed_save))
         if cmd.count_atoms(sele,quiet=1):
            cmd.edit(sele)
         cmd.delete(tmp_ed_save)
      elif cmd.count_atoms("((%s) and elem h)"%selection,quiet=1):
         print " Error: please pick a nitrogen or carbonyl carbon to grow from."
         raise QuietException
     
   cmd.delete(tmp_editor)
      


