# pymol -c generate2.py

# NOTE: obsolete -- PyMOL now uses Dunbrack rotamers by default

from chempy import io
from glob import glob
from copy import deepcopy

sys.path.append(".")

sc_raw = io.pkl.fromFile("sc_raw.pkl")

cutoff_angle = 35.0

sc_clus = {}

# this loop clusters and averages sidechain conformations

for resn in sc_raw.keys():
   resn_list = sc_raw[resn]
   resn_key_dict = {}
   # find the most commonly encountered set of torsions
   # (it is most likely to be correct)
   for set in resn_list:
      lst = deepcopy(set.keys())
      lst.sort()
      tup = tuple(lst)
      if resn_key_dict.has_key(tup):
         resn_key_dict[tup] = resn_key_dict[tup] + 1
      else:
         resn_key_dict[tup] = 1
   key_lst = []
   for a in resn_key_dict.keys():
      key_lst.append((resn_key_dict[a],a))
   key_lst.sort()
   resn_key = list(key_lst[-1][1])
   resn_key.sort()
   print resn,resn_key
   n_dihe = len(resn_key)
   print resn,len(sc_raw[resn])#,resn_key
   if n_dihe: # not glycine or alanine
      # list of dictionaries
      # [ {(...)=avg1, (...)=avg2, ... }, {(..)=avg1, ... }, ... ]
      avg_ang = []
      # list of list of dictionaries:
      # [ [ {(...)=ang1, (...)=ang2, ... }, {(..)=ang1, ... }, ... ] ]
      # where each row corresponds to the
      # averages used to generate the avg_ang [ [ {(..)=}, 
      all_ang = []
      while len(resn_list): # do we have any more cases to consider?
         # yes, lets consider the conformation "cur"
         cur = resn_list.pop()
         if len(cur)!=n_dihe:
#            print "skipping...",cur
            continue
         flag=1
         for k in resn_key:
            if not cur.has_key(k):
               flag=0
         if not flag:
            continue
         recomp_avg = None
         closest = None
         min_dev = 361.00
         cnt = 0
         for a in avg_ang: # we're going to compare it against all knowns
            max_dev = 0.0
            avg_dev = 0.0
            for k in resn_key: # for each known, compare all torsions
               dev = abs(cur[k]-a[k])
               if dev>180.0:
                  dev = 360.0 - dev
               avg_dev = avg_dev + dev
               if dev>max_dev:
                  max_dev = dev
            avg_dev = avg_dev / n_dihe
            if max_dev<cutoff_angle: # if this is acceptable...
               if min_dev>avg_dev: # and it is the closest, then remember
                  closest = cnt
                  min_dev = dev
            cnt = cnt + 1
         # we've now compared "cur" against all knowns...
         # so does it fit within an existing cluster?
         if closest!=None:
            # yes, it does, so add it into the group
            all_ang[closest].append(cur)
            recomp_avg = closest
         else:
            # nope, it is new
            all_ang.append([cur])
            avg_ang.append({})
            recomp_avg = len(avg_ang)-1
         # recomp an average conformation for the next round
         if recomp_avg!=None:
            avg_ang[recomp_avg]=deepcopy(all_ang[recomp_avg][0])
            for k in resn_key:
               avg_dev = 0.0
               for a in all_ang[recomp_avg]:
                  dev = cur[k]-a[k] 
                  if dev>180.0: dev = 360.0 - dev
                  if dev<-180.0: dev = 360.0 + dev
                  avg_dev = avg_dev + dev
               # compute average deviation, and save
               avg_deg = avg_dev / len(all_ang[recomp_avg])
               avg_ang[recomp_avg][k] = avg_ang[recomp_avg][k] + dev
      # now go in and compute frequency of occurrence, and sort
      tot = 0.0
      for a in all_ang:
         tot = tot + len(a)
      cnt = 0
      sort_ang = []
      for a in avg_ang:
         freq = len(all_ang[cnt])/tot
         a['FREQ'] = freq         
         sort_ang.append((freq,a))
         cnt = cnt + 1
      sort_ang.sort()
      avg_ang = []
      for a in sort_ang:
         avg_ang.insert(0,a[1])
      # at this point, we have a list of clustered, averaged torsions 
      # which can be used by the sidechain placement algorithm
      sc_clus[resn] = avg_ang
      print " reduced to:",len(sc_clus[resn])
      
io.pkl.toFile(sc_clus,"sc_library.pkl")


