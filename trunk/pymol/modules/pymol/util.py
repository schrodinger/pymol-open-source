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

import cmd
import math
import string
import pymol
from pymol import movie
# legacy mappings, remove in PyMOL 2.0

mload = movie.load
mrock = movie.rock
mroll = movie.roll

def cbag(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("carbon","(elem C and "+s+")")

def cbac(selection):
   s = str(selection)
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("cyan","(elem C and "+s+")")

def cbay(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("yellow","(elem C and "+s+")")

def cbas(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("salmon","(elem C and "+s+")")

def cbap(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("purple","(elem C and "+s+")")

def cbak(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("pink","(elem C and "+s+")")

def cbaw(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("hydrogen","(elem C and "+s+")")

def cbab(selection):
   s = str(selection)   
   cmd.color("magenta","("+s+")")
   cmd.color("oxygen","(elem O and "+s+")")
   cmd.color("nitrogen","(elem N and "+s+")")
   cmd.color("sulfer","(elem S and "+s+")")
   cmd.color("hydrogen","(elem H and "+s+")")
   cmd.color("slate","(elem C and "+s+")")

def performance(mode):
   mode = int(mode)
   if mode==0: # maximum quality
      cmd.set('line_smooth',1)
      cmd.set('depth_cue',1)         
      cmd.set('specular',1)
      cmd.set('surface_quality',1)
      cmd.set('stick_quality',15)
      cmd.set('sphere_quality',2)
      cmd.set('cartoon_sampling',14)
      cmd.set('ribbon_sampling',10)
      cmd.do("rebuild")
   elif mode==33:
      cmd.set('line_smooth',1)         
      cmd.set('depth_cue',1)         
      cmd.set('specular',1)
      cmd.set('surface_quality',0)
      cmd.set('stick_quality',8)
      cmd.set('sphere_quality',1)
      cmd.set('cartoon_sampling',7)
      cmd.do("rebuild")
   elif mode==66: # good perfomance
      cmd.set('line_smooth',0)
      cmd.set('depth_cue',0)         
      cmd.set('specular',1)
      cmd.set('surface_quality',0)
      cmd.set('stick_quality',8)
      cmd.set('sphere_quality',1)
      cmd.set('cartoon_sampling',6)
      cmd.do("rebuild")         
   else: # maximum performance
      cmd.set('line_smooth',0)
      cmd.set('depth_cue',0)
      cmd.set('specular',0)
      cmd.set('surface_quality',-1) # new
      cmd.set('stick_quality',5)
      cmd.set('sphere_quality',0)
      cmd.set('cartoon_sampling',3)
      cmd.do("rebuild")         

def hide_sele():
   arg = cmd.get_names("selections")
   for a in arg:
      cmd.disable(a)

def hbond(a,b,cutoff=3.3):
   st = "(%s and (%s around %4.2f) and elem N,O),(%s and (%s around %4.2f) and elem N,O),%4.2f" % (a,b,cutoff,b,a,cutoff,cutoff)
#   cmd.dist("hbond",st)

def cbc(selection='(all)'): # NOT THREAD SAFE
   '''
   Color all chains a different color
   '''
   pymol.stored.chain = {}
   cmd.iterate("(%s)"%selection,"stored.chain[chain]=1")
   c = 7
   for a in pymol.stored.chain.keys():
      if len(a):
         print ("%d,(chain %s)"%(c,a))
         cmd.color("%d"%c,"(chain %s)"%a)
         c = c + 1

color_chains = cbc

def sum_charge(*arg): # NOT THREAD SAFE
   result = None
   try:
      obj = "all"
      if len(arg):
         obj = arg

      pymol.stored._sum_charge = 0.0
      cmd.iterate("(%s)"%obj,
                  "stored._sum_charge=stored._sum_charge+partial_charge")
      result = pymol.stored._sum_charge
      print " sum_charge: %6.4f"%result
   except:
      print " sum_charge: an error occurred."
   return result

def ray_shadows(mode):
   if mode=='light': # maximum quality
      cmd.set('power',3.0)
      cmd.set('ambient',0.3)
      cmd.set('spec_power',40)
      cmd.set('reflect',1.2)         
      cmd.set('direct',0.35)
      cmd.set('gamma',1.2)
   elif mode=='matte':
      cmd.set('power',1.5) 
      cmd.set('spec_power',20)
      cmd.set('ambient',0.12)
      cmd.set('reflect',0.9) 
      cmd.set('direct',0.25)
      cmd.set('gamma',1.30)
   elif mode=='medium':
      cmd.set('power',1.0) # 0.7
      cmd.set('spec_power',60) # was 50
      cmd.set('ambient',0.12)
      cmd.set('reflect',0.9) 
      cmd.set('direct',0.25) 
      cmd.set('gamma',1.3) 
   elif mode=='heavy':
      cmd.set('power',0.3) 
      cmd.set('spec_power',90) # was 60
      cmd.set('ambient',0.08)
      cmd.set('reflect',0.65) # was 0.75         
      cmd.set('direct',0.06)
      cmd.set('gamma',1.4) # was 1.5
   
def ff_copy(src,dst): # NOT THREAD SAFE
   pymol._rcopy = pymol.Scratch_Storage()
   pymol._rcopy.pc={}
   pymol._rcopy.tt={}
   cmd.iterate("(%s)"%src,"_rcopy.pc[name]=partial_charge")
   cmd.alter("(%s)"%dst,"partial_charge=_rcopy.pc[name]")
   cmd.iterate("(%s)"%src,"_rcopy.tt[name]=text_type")
   cmd.alter("(%s)"%dst,"text_type=_rcopy.tt[name]")
   del pymol._rcopy
   
def b2vdw(*arg):
   if not len(arg):
      sele = 'all'
   else:
      sele = arg[0]
   # use B values to create RMS VDW spheres
   # rms = sqrt(b/(8*(PI^2)))
   cmd.alter("(%s)"%sele,"vdw=math.sqrt(b/78.9568352087)")
   
def phipsi(selection="(pk1)"): # NOT THREAD SAFE
   n_sele =   "((byres (%s)) & name n)"%selection
   c_sele =   "((byres (%s)) & name c)"%selection
   ca_sele =  "((byres (%s)) & name ca)"%selection
   cm_sele = "((neighbor (%s)) and not (byres (%s)))"%(n_sele,n_sele)
   np_sele = "((neighbor (%s)) and not (byres (%s)))"%(c_sele,c_sele)
   cmd.feedback("push")
   cmd.feedback("disable","selector","everythin")
   cm_cnt = cmd.select("_pp_cm",cm_sele)
   n_cnt = cmd.select("_pp_n",n_sele)
   c_cnt = cmd.select("_pp_c",c_sele)
   ca_cnt = cmd.select("_pp_ca",ca_sele)
   np_cnt = cmd.select("_pp_np",np_sele)
   if(cm_cnt and n_cnt and ca_cnt and c_cnt):
      phi = cmd.get_dihedral("_pp_c","_pp_ca","_pp_n","_pp_cm")
   else:
      phi = None
   if(n_cnt and ca_cnt and c_cnt and np_cnt):
      psi = cmd.get_dihedral("_pp_np","_pp_c","_pp_ca","_pp_n")
   else:
      psi = None
   cmd.feedback("pop")
   cmd.delete("_pp_cm")
   cmd.delete("_pp_n")
   cmd.delete("_pp_c")
   cmd.delete("_pp_ca")
   cmd.delete("_pp_np")
   return (phi,psi)

def rainbow(selection="(name ca and alt '',A)",reverse=0): # NOT THREAD SAFE

   cmd.feedback("push")
   cmd.feedback("disable","executive","actions")

   # your basic rainbow...
   
   list = [
      (0,0,255),      
      (0,0,255),
      (0,128,255),
      (0,255,255),
      (0,255,128),           
      (0,255,0),
      (128,255,0),
      (255,255,0),
      (255,128,0),
      (255,0,0),
      (255,0,0)      
      ]
   if reverse:
      list.reverse()
   #
   last = list.pop(0)
   cmd.set_color("_000",[last[0]/255.0,last[1]/255.0,last[2]/255.0])
   c = 1
   for a in list:
      for b in range(1,21):
         b0 = b/20.0
         b1 = 1.0-b0
         cname = "_%03d"%c
         r = last[0]*b1+a[0]*b0
         g = last[1]*b1+a[1]*b0
         b = last[2]*b1+a[2]*b0
         cmd.set_color(cname,[r/255.0,g/255.0,b/255.0])
         c = c + 1
      last = a

   cas = cmd.index("((byres ("+selection+")) and name ca and not het)")
   l = len(cas)
   if not len(cas):
      return
   c = 0
   for a in cas:
      col = int((200*c)/l)
      cmd.color("_%03d"%col,"((%s) and (byres %s`%d))"%(selection,a[0],a[1]))
      c = c + 1

   cmd.feedback("pop")
   
def ss(selection="(name ca and alt '',A)",state=1): # NOT THREAD SAFE

   print ' util.ss: WARNING: This is not a "correct" secondary structure'
   print ' util.ss: assignment algorithm!  Please use only as a last resort.'
   
   cmd.feedback("push")
   cmd.feedback("disable","executive","actions")
   
   ss_pref = "_sss"
   sss1 = ss_pref+"1"
   cnt = cmd.select(sss1,"((byres ("+selection+")) and name ca and not het)")
   print " util.ss: initiating secondary structure assignment on %d residues."%cnt
   cas = cmd.index(sss1)
   if not len(cas):
      return
   # set cartoon mode to auto over the selection
   
   cmd.cartoon("auto",sss1)

   print " util.ss: extracting sequence and relationships..."

   # get CA list
   
   res_list = []
   pymol._ss = pymol.Scratch_Storage()
   pymol._ss.res_list = res_list
   cmd.iterate(sss1,'_ss.res_list.append((model,index))')

   # generate atom-to-residue conversion dictionaries

   ca_dict = {}
   n_dict = {}
   o_dict = {}
   scr_dict = {} # scr = segment,chain,resi 
   pymol._ss.n_dict = n_dict
   pymol._ss.o_dict = o_dict
   pymol._ss.scr_dict = scr_dict
   pymol._ss.ca_dict = ca_dict
   cmd.iterate(sss1,
               '_ss.scr_dict[(model,index)]=(segi,chain,resi)') # CA's
   cmd.iterate("((byres "+sss1+") and n;n)"
               ,'_ss.scr_dict[(model,index)]=(segi,chain,resi)') # N's
   cmd.iterate("((byres "+sss1+") and n;o)",
               '_ss.scr_dict[(model,index)]=(segi,chain,resi)') # O's
   cmd.iterate(sss1,
               '_ss.ca_dict[(segi,chain,resi)] = (model,index)')
   cmd.iterate("((byres "+sss1+") and n;n)",
               '_ss.n_dict[(segi,chain,resi)] = (model,index)')
   cmd.iterate("((byres "+sss1+") and n;o)",
               '_ss.o_dict[(segi,chain,resi)] = (model,index)')

   scr_dict[None]=None
   o_dict[None]=None
   n_dict[None]=None
   ca_dict[None]=None
   
   # create special version of cas with gaps

   gap = [None,None,None,None]  
   # gap large enough to distinguish i+4 interations from gaps
   last = None
   for a in res_list:
      if last!=None:
         if(cmd.count_atoms(
            "((neighbor(neighbor(neighbor (%s`%d)))) and (%s`%d))"%
            (last[0],last[1],a[0],a[1]),quiet=1)==0):
            gap.extend([None,None,None,None])
      gap.append(a)
      last = a
   gap.extend([None,None,None,None])

   print " util.ss: analyzing phi/psi angles (slow)..."

   # generate reverse-lookup for gap indices

   ss = {}

   c = 0
   gap_idx = {}
   for a in gap:
      gap_idx[a] = c
      c = c + 1

   # secondary structure database...
   
   ss = {}
   ss[None]=None
   
   # make decisions based on phi/psi

   for a in cas:
      ss[a] = 'L' # default
   phipsi = cmd.get_phipsi(sss1,state)
   for a in phipsi.keys():
      (phi,psi) = phipsi[a]
#      print scr_dict[a],(phi,psi)
      if (phi!=None) and (psi!=None):
         if ((phi<-45) and (phi>-160) and
             (psi<-170) or (psi>10)): # beta?
            ss[a] = 's'
         elif ((phi<-45) and (phi>-160) and
               (psi>-80) and (psi<-25)): # helix?
            ss[a] = 'H'
            
   print " util.ss: finding hydrogen bonds..."
   
   # find all pairwise hydrogen bonds and make note of them in dict

   hb = cmd.find_pairs("((byres "+sss1+") and n;n)",
                       "((byres "+sss1+") and n;o)",mode=1,
                       cutoff=3.7,angle=55,
                       state1=state,state2=state)
   
   hb_dict = {}  # [((N-atom) (O-atom))] = 1
   n_hb_dict = {} # [(N-atom)] = [(O-atom),...]
   o_hb_dict = {} # [(O-atom)] = [(N-atom),...]
   for a in hb:
#      cmd.dist("(%s`%d)"%a[0],"(%s`%d)"%a[1])
      hb_dict[a] = 1
      n = a[0]
      o = a[1]
      if not n_hb_dict.has_key(n): n_hb_dict[n]=[]
      if not o_hb_dict.has_key(o): o_hb_dict[o]=[]
      n_hb_dict[n].append(o)
      o_hb_dict[o].append(n)

   # check to insure that all helical residues have at least an i +/- 4
   # hydrogen bond

   for c in xrange(4,len(gap)-4):
      a = gap[c]
      if ss[a]=='H':
         aN = n_dict[scr_dict[a]]
         aO = o_dict[scr_dict[a]]
         am4O = o_dict[scr_dict[gap[c-4]]]
         ap4N = n_dict[scr_dict[gap[c+4]]]
         if not hb_dict.has_key((aN,am4O)):
            if not hb_dict.has_key((ap4N,aO)):
               ss[a]='L'

   print " util.ss: verifying beta sheets..."
   
   # check to insure that all beta residues have proper interactions

   rep_dict = {}
   repeat = 1
   while repeat:
      repeat = 0
      c = 4
      cc = len(gap)-4
      while c<cc:
         a1 = gap[c]
         if (ss[a1] in ['s','S']) and not rep_dict.has_key(a1):
            rep_dict[a1] = 1
            valid = 0
            scr_a1 = scr_dict[a1]
            # look for antiparallel 2:2 H-bonds (NH-O=C + C=O-HN) 
            n_a1_atom = n_dict[scr_a1]
            o_a1_atom = o_dict[scr_a1]
            if (n_hb_dict.has_key(n_a1_atom) and 
                o_hb_dict.has_key(o_a1_atom)):
               for n_hb_atom in n_hb_dict[n_a1_atom]:
                  for o_hb_atom in o_hb_dict[o_a1_atom]:
                     n_hb_scr = scr_dict[n_hb_atom]
                     o_hb_scr = scr_dict[o_hb_atom]
                     if o_hb_scr == n_hb_scr:
                        b1 = ca_dict[o_hb_scr]
                        if abs(c-gap_idx[b1])>2:
                           ss[b1] = 'S' 
                           ss[a1] = 'S' 
                           valid = 1
            # look for antiparallel offset HB (i,i+2,j,j-2)
            a3 = gap[c+2]
            if (a3!=None):
               scr_a3 = scr_dict[a3]
               o_a1_atom = o_dict[scr_a1]
               n_a3_atom = n_dict[scr_a3]
               if (n_hb_dict.has_key(n_a3_atom) and
                   o_hb_dict.has_key(o_a1_atom)):               
                  for n_hb_atom in n_hb_dict[n_a3_atom]:
                     for o_hb_atom in o_hb_dict[o_a1_atom]:
                        n_hb_scr = scr_dict[n_hb_atom]
                        o_hb_scr = scr_dict[o_hb_atom]
                        b1 = ca_dict[o_hb_scr]
                        if b1!=None:
                           b1_i = gap_idx[b1]
                           if abs(c-b1_i)>2: # no turns!
                              b3 = gap[b1_i-2]
                              if b3!=None:
                                 b3_scr = scr_dict[b3]
                                 if b3_scr == n_hb_scr:
                                    a2 = gap[c+1]
                                    b2 = gap[gap_idx[b1]-1]
                                    ss[b1] = 'S'
                                    ss[b3] = 'S'
                                    ss[a1] = 'S'
                                    ss[a3] = 'S'
                                    if ss[a2]=='L': ss[a2] = 's'
                                    if ss[b2]=='L': ss[b2] = 's'
                                    valid = 1
            # look for antiparallel offset HB (i,i-2,j,j+2)
            a3 = gap[c-2]
            if (a3!=None):
               scr_a3 = scr_dict[a3]
               n_a1_atom = n_dict[scr_a1]
               o_a3_atom = o_dict[scr_a3]
               if (n_hb_dict.has_key(n_a1_atom) and
                   o_hb_dict.has_key(o_a3_atom)):               
                  for n_hb_atom in n_hb_dict[n_a1_atom]:
                     for o_hb_atom in o_hb_dict[o_a3_atom]:
                        n_hb_scr = scr_dict[n_hb_atom]
                        o_hb_scr = scr_dict[o_hb_atom]
                        b1 = ca_dict[o_hb_scr]
                        if b1!=None:
                           b1_i = gap_idx[b1]
                           if abs(c-b1_i)>2: # no turns!
                              b3 = gap[b1_i-2]
                              if b3!=None:
                                 b3_scr = scr_dict[b3]
                                 if b3_scr == n_hb_scr:
                                    a2 = gap[c-1]
                                    b2 = gap[gap_idx[b1]-1]
                                    ss[b1] = 'S'
                                    ss[b3] = 'S'
                                    ss[a1] = 'S'
                                    ss[a3] = 'S'
                                    if ss[a2]=='L': ss[a2] = 's'
                                    if ss[b2]=='L': ss[b2] = 's'
                                    valid = 1
            # look for parallel 1:3 HB (i,j-1,j+1)
            n_a1_atom = n_dict[scr_a1]
            o_a1_atom = o_dict[scr_a1]
            if (n_hb_dict.has_key(n_a1_atom) and
                o_hb_dict.has_key(o_a1_atom)):
               for n_hb_atom in n_hb_dict[n_a1_atom]:
                  for o_hb_atom in o_hb_dict[o_a1_atom]:
                     n_hb_scr = scr_dict[n_hb_atom]
                     o_hb_scr = scr_dict[o_hb_atom]
                     b0 = ca_dict[n_hb_scr]
                     if b0!=None:
                        b2 = gap[gap_idx[b0]+2]
                        if b2!=None:
                           b2_scr = scr_dict[b2]
                           if b2_scr == o_hb_scr:
                              b1 = gap[gap_idx[b0]+1]
                              ss[a1] = 'S' 
                              ss[b0] = 'S'
                              if ss[b1]=='L': ss[b1]='s'
                              ss[b2] = 'S'
                              valid = 1
                              repeat = 1
            if not valid:
               ss[a1] = 'L'
         c = c + 1

   # automatically fill 1 residue gaps in helices and well-defined sheets
   c = 4
   cc = len(gap)-6
   while c<cc:
      a1 = gap[c]
      a3 = gap[c+2]
      ss_a1 = ss[a1]
      ss_a3 = ss[a3]
      if (ss_a1==ss_a3) and (ss_a1 in ['S','H']):
         a2 = gap[c+1]
         ss[a2] = ss_a1
      c = c + 1

   # remove singleton sheet residues
   c = 4
   cc = len(gap)-4
   while c<cc:
      a0 = gap[c-1]
      a1 = gap[c]
      a2 = gap[c+1]
      if ss[a1] in ['s','S']:
         if ((not ss[a0] in ['s','S']) and
             (not ss[a2] in ['s','S'])):
             ss[a1] = 'L'
      c = c + 1

   # remove sheet residues which aren't next to another sheet 
   c = 4
   cc = len(gap)-4
   while c<cc:
      a1 = gap[c]
      if ss[a1]=='S':
         a1 = gap[c]
         scr_a1 = scr_dict[a1]
         # look for hydrogen bonds to another sheet
         n_a1_atom = n_dict[scr_a1]
         o_a1_atom = o_dict[scr_a1]
         certain = 0
         if n_hb_dict.has_key(n_a1_atom):
            for n_hb_atom in n_hb_dict[n_a1_atom]:
               n_hb_ca_atom=ca_dict[scr_dict[n_hb_atom]]
               if ss[n_hb_ca_atom]=='S':
                  certain = 1
                  break
         if o_hb_dict.has_key(o_a1_atom):
            for o_hb_atom in o_hb_dict[o_a1_atom]:
               o_hb_ca_atom=ca_dict[scr_dict[o_hb_atom]]
               if ss[o_hb_ca_atom]=='S':
                  certain = 1
                  break
         if not certain:
            ss[a1] = 's'
      c = c + 1

   # remove questionable sheet residues
   c = 4
   cc = len(gap)-4
   while c<cc:
      a0 = gap[c-1]
      a1 = gap[c]
      a2 = gap[c+1]
      if ss[a1]=='s':
         if (not ((ss[a0]=='S') and (ss[a2]=='S'))):
            ss[a1] = 'L'
      c = c + 1

   # extend helices if hydrogen bonding requirements are met
   rep_dict = {}
   repeat = 1
   while repeat:
      repeat = 0
      c = 4
      cc = len(gap)-4
      while c<cc:
         a = gap[c]
         if not rep_dict.has_key(a):
            if ss[gap[c+1]]=='H':
               rep_dict[a] = 1
               if ss[a]!='H': # N-terminal end
                  aO = o_dict[scr_dict[a]]
                  ap4N = n_dict[scr_dict[gap[c+4]]]
                  ap3N = n_dict[scr_dict[gap[c+3]]]
                  if hb_dict.has_key((ap4N,aO)) or hb_dict.has_key((ap3N,aO)):
                     ss[a]='H'
                     repeat = 1
                     c = c - 5
                     if c<4: c=4
            if ss[gap[c-1]]=='H':
               a = gap[c]
               if ss[a]!='H': # C-terminal end
                  rep_dict[a] = 1
                  aN = n_dict[scr_dict[a]]
                  am4O = o_dict[scr_dict[gap[c-4]]]
                  am3O = o_dict[scr_dict[gap[c-3]]]
                  if hb_dict.has_key((aN,am4O)) or hb_dict.has_key((aN,am3O)):
                     ss[a]='H'
                     repeat = 1
                     c = c - 5
                     if c<4: c=4
         c = c + 1

   # remove doubleton helices

   c = 4
   cc = len(gap)-5
   while c<cc:
      a0 = gap[c-1]
      a1 = gap[c]
      a2 = gap[c+1]
      a3 = gap[c+2]
      ss_a0 = ss[gap[c-1]]
      ss_a1 = ss[gap[c]]
      ss_a2 = ss[gap[c+1]]
      ss_a3 = ss[gap[c+2]]
      if ss_a1=='H':
         if (ss_a2==ss_a1) and (ss_a0!=ss_a2) and (ss_a2!=ss_a3):
            ss[a1] = 'L'
            ss[a2] = 'L'
      c = c + 1

   # remove totally unreasonable helix and sheet residues

   c = 4
   cc = len(gap)-5
   while c<cc:
      a1 = gap[c]
      ss_a1 = ss[gap[c]]
      if ss_a1=='H':
         if phipsi.has_key(a1):
            (phi,psi) = phipsi[a1]
            if (phi>0) and (phi<150):
               ss[a1] = 'L'
            elif((psi<-120) or (psi>140)):
               ss[a1] = 'L'
      elif ss_a1 in ['S','s']:
         if phipsi.has_key(a1):
            (phi,psi) = phipsi[a1]
            if (phi>45) and (phi<160):
               ss[a1] = 'L'
#            if (psi<-30) and (psi>-150):
            if (psi<-65) and (psi>-150):
               ss[a1] = 'L'
         
      c = c + 1


   for x in range(1,3):
      # remove singleton sheet residues
      c = 4
      cc = len(gap)-4
      while c<cc:
         a0 = gap[c-1]
         a1 = gap[c]
         a2 = gap[c+1]
         if ss[a1] in ['s','S']:
            if ((not ss[a0] in ['s','S']) and
                (not ss[a2] in ['s','S'])):
                ss[a1] = 'L'
         c = c + 1

      # remove sheet residues which aren't next to another sheet 
      c = 4
      cc = len(gap)-4
      while c<cc:
         a1 = gap[c]
         if ss[a1]=='S':
            a1 = gap[c]
            scr_a1 = scr_dict[a1]
            # look for hydrogen bonds to another sheet
            n_a1_atom = n_dict[scr_a1]
            o_a1_atom = o_dict[scr_a1]
            certain = 0
            if n_hb_dict.has_key(n_a1_atom):
               for n_hb_atom in n_hb_dict[n_a1_atom]:
                  n_hb_ca_atom=ca_dict[scr_dict[n_hb_atom]]
                  if ss[n_hb_ca_atom]=='S':
                     certain = 1
                     break
            if o_hb_dict.has_key(o_a1_atom):
               for o_hb_atom in o_hb_dict[o_a1_atom]:
                  o_hb_ca_atom=ca_dict[scr_dict[o_hb_atom]]
                  if ss[o_hb_ca_atom]=='S':
                     certain = 1
                     break
            if not certain:
               ss[a1] = 's'
         c = c + 1

      # remove questionable sheet residues
      c = 4
      cc = len(gap)-4
      while c<cc:
         a0 = gap[c-1]
         a1 = gap[c]
         a2 = gap[c+1]
         if ss[a1]=='s':
            if (not ((ss[a0]=='S') and (ss[a2]=='S'))):
               ss[a1] = 'L'
         c = c + 1

#      lst = ss.keys()
#      lst.sort()
#      for a in lst: print scr_dict[a],ss[a]
      
   # assign protein
   for a in cas:
      if ss[a]=='s':
         ss[a]='S'
      
   cmd.alter(sss1,"ss ='L'")
   for a in cas:
      if ss[a]!='L':
         cmd.alter("(%s`%d)"%a,"ss='%s'"%ss[a])

   cmd.feedback("pop")

   del pymol._ss # IMPORTANT
   cmd.delete(sss1)
   
   #
#   print conn_hash.keys()
   print " util.ss: assignment complete."

