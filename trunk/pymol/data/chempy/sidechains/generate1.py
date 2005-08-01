# pymol -c generate1.py

# NOTE: obsolete -- PyMOL now uses Dunbrack rotamers by default

from chempy import io
from glob import glob

sys.path.append(".")

from pymol import cmd

import pymol
from pymol import cmd

def find_torsions(res_sele):
   # this routine gets a list of atom ids which correspond to the
   # free torsions in the sidechain
   res_sele = "("+res_sele+")"
   if cmd.count_atoms("(%s and n;N,CA)"%res_sele,quiet=1)<2:
      return []
   result = []
   done = 0
   # set up the starting trunk
   cmd.select("trunk","(%s and n;N,CA,C,O)"%res_sele)
   first_id = cmd.identify("(%s and n;N)"%res_sele)[0]
   second_id = cmd.identify("(%s and n;CA)"%res_sele)[0]
   to_do_list = [(first_id,second_id)]
   while len(to_do_list):
      (first_id,second_id) = to_do_list.pop(0)
      # find candidates for the 3rd atom in the dihedral
      cmd.select("candi","(%s and (neighbor trunk))"%res_sele)
      pymol.stored.list = []
      cmd.iterate("candi","stored.list.append((name,ID))")
      candi_list = pymol.stored.list
      if len(candi_list):
         # we may have a dihedral
         dihe_list = []
         for a in candi_list:
            if cmd.select("termi","(%s and (not trunk) and (neighbor id %d))"%
                          (res_sele,a[1]))>0:
               # apparently we do, but is it acyclic
               # (i.e., can edit split the molecule?)
               cmd.edit("(%s and trunk and (neighbor id %d))"%(res_sele,a[1]),
                        "(%s and id %d)"%(res_sele,a[1]))
               if 'pkfrag2' in cmd.get_names('selections'):
                  # yes, so this is a valid dihedral
                  dihe_list.append(a)
               else:
                  # no, so add the third atom to the trunk
                  cmd.select("trunk","(trunk or (%s and id %d))"%(res_sele,a[1]))                  
            else:
               # no fourth atom, so add this atom to the trunk
               cmd.select("trunk","(trunk or (%s and id %d))"%(res_sele,a[1]))
         if len(dihe_list):
            # choose 3rd atom using alphanumeric order
            dihe_list.sort()
            third_id = dihe_list[0][1]
            # if there is another third atom, then repeat later      
            if len(dihe_list)>1: 
               to_do_list.insert(0,(first_id,second_id))
            # now choose the 4th atom, which we know exists, using a similar criterion
            cmd.select("termi","(%s and (not trunk) and (neighbor id %d))"%
                       (res_sele,third_id))
            pymol.stored.list=[]
            cmd.iterate("termi","stored.list.append((name,ID))")
            termi_list = pymol.stored.list
            termi_list.sort()
            fourth_id = termi_list[0][1]
            # at this point, we should have a complete dihedral
            # add the third atom into the trunk, and store the second and third for
            # outward extension later on
            to_do_list.append((second_id,third_id))
            cmd.select("trunk","(trunk or (%s and id %d))"%(res_sele,third_id))      
#            cmd.show('sticks','(%s and (id %d,%d,%d,%d))'%
#                     (res_sele,first_id,second_id,third_id,fourth_id))
            result.append((first_id,second_id,third_id,fourth_id))
   return result

def measure_all_torsions(sele,prot_dict):
   id_list = cmd.identify("((%s) and name ca)"%sele)
   for a in id_list:
      res_sele = "(alt '' and (byres id %d) and (%s))"%(a,sele)
      lst = find_torsions(res_sele)
      if len(lst):
         cmd.iterate("(%s and n;ca)"%res_sele,"stored.resn=resn")
         resn = pymol.stored.resn
         print resn
         name_list = []
         for a in lst:
            pymol.stored.names = []
            for b in a:
               cmd.iterate("(%s and id %d)"%(res_sele,b),"stored.names.append(name)")
            name_list.append(tuple(pymol.stored.names))
         res_dict = {}
         for a in name_list:
            dihe = cmd.get_dihedral("(%s and name %s)"%(res_sele,a[0]),
                                    "(%s and name %s)"%(res_sele,a[1]),
                                    "(%s and name %s)"%(res_sele,a[2]),
                                    "(%s and name %s)"%(res_sele,a[3]))
            res_dict[a]=dihe
         if not prot_dict.has_key(resn):
            prot_dict[resn]=[]
         prot_dict[resn].append(res_dict)

sc_dict = {}

cmd.feedback('disable','selector','result')
cmd.feedback('disable','executive','action')

for fyle in glob("source/*"):
   cmd.load(fyle)
   cmd.remove("(het)")
   cmd.remove("(hydro or not alt '')")
   cmd.remove("(byres b>25)") # only use well-defined residues
   cmd.remove("(not (byres name ca))")
   cmd.remove("(not (byres name n))")
   cmd.remove("(not (byres name c))")   
   cmd.refresh()
   measure_all_torsions("all",sc_dict)
   cmd.delete("all")
   
#for a in sc_dict.keys():
#   print a,dict[a]
#print len(dict.keys())
      
io.pkl.toFile(sc_dict,"sc_raw.pkl")

